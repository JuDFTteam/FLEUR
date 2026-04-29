#!/usr/bin/env python3
"""Discover SLURM machine/partition properties and build MachineConfig JSON.

This module queries SLURM via ``sinfo`` and ``scontrol`` and creates a
partition-centric machine configuration compatible with ``MachineConfig``.

Only ``name`` and ``description`` are machine-level fields. All scheduler and
resource settings are stored per partition.
"""

from __future__ import annotations

import json
import math
import re
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .slurm_generator import MachineConfig, Partition, get_machines_directory


DEFAULT_COMMAND = "srun fleur_MPI"


@dataclass
class PartitionDiscovery:
    """Intermediate representation of discovered partition settings."""

    name: str
    default: bool = False
    account: Optional[str] = None
    max_nodes: int = 1
    max_runtime: str = "24:00:00"
    cores_per_node: int = 1
    gpus_per_node: int = 0
    memory_per_node_gb: int = 1
    gpu_type: Optional[str] = None
    gpu_syntax: str = "gpus"
    use_mem_option: bool = True
    modules_needed: list[str] = None
    command: Optional[str] = None
    shell_commands: list[str] = None
    description: str = ""

    def to_partition(self) -> Partition:
        """Convert to Partition dataclass used by MachineConfig."""
        return Partition(
            name=self.name,
            default=self.default,
            account=self.account,
            max_nodes=self.max_nodes,
            max_runtime=self.max_runtime,
            cores_per_node=self.cores_per_node,
            gpus_per_node=self.gpus_per_node,
            memory_per_node_gb=self.memory_per_node_gb,
            gpu_type=self.gpu_type,
            gpu_syntax=self.gpu_syntax,
            use_mem_option=self.use_mem_option,
            modules_needed=list(self.modules_needed or []),
            command=self.command,
            shell_commands=list(self.shell_commands or []),
            description=self.description,
        )


class SlurmDiscoveryError(RuntimeError):
    """Raised when SLURM discovery fails."""


def _run_command(args: list[str], timeout: int = 20) -> str:
    """Run command and return stdout; raise with stderr on failure."""
    try:
        result = subprocess.run(
            args,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        cmd = shlex.join(args)
        raise SlurmDiscoveryError(
            f"Command not found: {cmd}\n"
            "SLURM discovery requires access to SLURM client tools (e.g. sinfo, scontrol).\n"
            "Run this function on the machine where the SLURM installation is accessible."
        ) from exc
    except subprocess.CalledProcessError as exc:
        cmd = shlex.join(args)
        stderr = (exc.stderr or "").strip()
        raise SlurmDiscoveryError(f"Command failed: {cmd}\n{stderr}") from exc
    except subprocess.TimeoutExpired as exc:
        cmd = shlex.join(args)
        raise SlurmDiscoveryError(f"Command timed out: {cmd}") from exc

    return result.stdout.strip()


def _parse_key_value_line(line: str) -> dict[str, str]:
    """Parse an scontrol -o style key=value line into a dict."""
    parsed: dict[str, str] = {}
    for token in line.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        parsed[key] = value
    return parsed


def _normalize_partition_name(name: str) -> tuple[str, bool]:
    """Return partition name and default flag from sinfo '%P' field."""
    clean = name.strip()
    is_default = clean.endswith("*")
    clean = clean.rstrip("*")
    return clean, is_default


def _parse_int_field(value: Optional[str], fallback: int) -> int:
    """Parse integer-like SLURM field with fallback for unlimited/unknown."""
    if not value:
        return fallback
    up = value.upper()
    if up in {"UNLIMITED", "INFINITE", "N/A", "NONE"}:
        return fallback

    # Handles values such as "128", "4+", "1-32" by taking the last integer.
    matches = re.findall(r"\d+", value)
    if not matches:
        return fallback
    return int(matches[-1])


def _parse_memory_mb_to_gb(value: Optional[str], fallback_gb: int = 1) -> int:
    """Convert sinfo '%m' memory (MB) to GB, rounded up."""
    if not value:
        return fallback_gb
    matches = re.findall(r"\d+", value)
    if not matches:
        return fallback_gb
    mb = int(matches[0])
    return max(1, int(math.ceil(mb / 1024.0)))


def _parse_gres_gpu(gres: Optional[str]) -> tuple[int, Optional[str]]:
    """Extract GPU count and GPU type from GRES string.

    Examples:
    - "gpu:4" -> (4, None)
    - "gpu:a100:4(S:0-3)" -> (4, "a100")
    - "gpu:v100:2" -> (2, "v100")
    """
    if not gres:
        return 0, None

    token = None
    for piece in gres.split(","):
        piece = piece.strip()
        if piece.startswith("gpu:"):
            token = piece
            break
    if token is None:
        return 0, None

    # Remove optional suffix in parentheses, e.g. "(S:0-3)"
    token = token.split("(", 1)[0]

    parts = token.split(":")
    if len(parts) == 2:
        # gpu:4
        try:
            return int(parts[1]), None
        except ValueError:
            return 0, None

    if len(parts) >= 3:
        # gpu:type:count or gpu:type:model:count
        try:
            count = int(parts[-1])
        except ValueError:
            count = 0
        gpu_type = ":".join(parts[1:-1]) if len(parts) > 2 else None
        gpu_type = gpu_type or None
        return count, gpu_type

    return 0, None


def _detect_cluster_name() -> str:
    """Detect cluster name from scontrol config."""
    out = _run_command(["scontrol", "show", "config", "-o"])
    cfg = _parse_key_value_line(out)
    return cfg.get("ClusterName", "slurm_cluster")


def _discover_partition_data() -> dict[str, PartitionDiscovery]:
    """Discover partition data from sinfo + scontrol.

    Returns:
        Mapping partition_name -> PartitionDiscovery
    """
    partitions: dict[str, PartitionDiscovery] = {}

    # sinfo gives concise partition-level resources from live scheduler state.
    # %P: partition name (* marks default), %c CPUs/node, %m memory MB/node,
    # %G GRES, %D total nodes in partition.
    sinfo_out = _run_command(["sinfo", "-h", "-o", "%P|%c|%m|%G|%D"])
    for raw_line in sinfo_out.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split("|", 4)
        if len(fields) != 5:
            continue

        p_name_raw, cpus_s, mem_s, gres_s, nodes_s = [f.strip() for f in fields]
        p_name, default_from_sinfo = _normalize_partition_name(p_name_raw)
        if not p_name:
            continue

        cpus = _parse_int_field(cpus_s, fallback=1)
        mem_gb = _parse_memory_mb_to_gb(mem_s, fallback_gb=1)
        gpus, gpu_type = _parse_gres_gpu(gres_s)
        nodes = _parse_int_field(nodes_s, fallback=1)

        existing = partitions.get(p_name)
        if existing is None:
            existing = PartitionDiscovery(name=p_name)
            partitions[p_name] = existing

        # When sinfo yields multiple lines (states), keep conservative maxima.
        existing.default = existing.default or default_from_sinfo
        existing.cores_per_node = max(existing.cores_per_node, cpus)
        existing.memory_per_node_gb = max(existing.memory_per_node_gb, mem_gb)
        existing.gpus_per_node = max(existing.gpus_per_node, gpus)
        existing.max_nodes = max(existing.max_nodes, nodes)
        if existing.gpu_type is None and gpu_type:
            existing.gpu_type = gpu_type

    # scontrol adds static limits like MaxNodes and MaxTime.
    scontrol_out = _run_command(["scontrol", "show", "partition", "-o"])
    for raw_line in scontrol_out.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        info = _parse_key_value_line(line)
        p_name = info.get("PartitionName", "").strip()
        if not p_name:
            continue

        p = partitions.get(p_name)
        if p is None:
            p = PartitionDiscovery(name=p_name)
            partitions[p_name] = p

        default_flag = info.get("Default", "NO").upper() in {"YES", "TRUE"}
        p.default = p.default or default_flag

        p.max_nodes = _parse_int_field(info.get("MaxNodes"), fallback=p.max_nodes)
        p.max_runtime = info.get("MaxTime", p.max_runtime)

        allow_accounts = info.get("AllowAccounts")
        if allow_accounts and allow_accounts.upper() not in {"ALL", "(NULL)", "NONE"}:
            # Pick first account from list if available.
            p.account = allow_accounts.split(",", 1)[0]

    if not partitions:
        raise SlurmDiscoveryError("No partitions discovered from sinfo/scontrol")

    # Ensure a deterministic default partition exists.
    if not any(p.default for p in partitions.values()):
        first_name = sorted(partitions.keys())[0]
        partitions[first_name].default = True

    return partitions


def discover_machine_config(
    machine_name: Optional[str] = None,
    description: Optional[str] = None,
    default_command: str = DEFAULT_COMMAND,
    modules_needed: Optional[list[str]] = None,
    shell_commands: Optional[list[str]] = None,
    gpu_syntax: str = "gres",
    use_mem_option: bool = True,
) -> MachineConfig:
    """Discover SLURM machine configuration and return a MachineConfig object.

    Args:
        machine_name: Optional machine name; if omitted derived from ClusterName.
        description: Optional machine description.
        default_command: Command template assigned to each partition.
        modules_needed: Optional modules assigned to each partition.
        shell_commands: Optional shell commands assigned to each partition.
        gpu_syntax: GPU directive syntax ("gpus" or "gres").
        use_mem_option: Whether generated configs should emit ``--mem``.
    """
    cluster_name = machine_name or _detect_cluster_name()
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    desc = description or f"Auto-discovered from SLURM on {timestamp}"

    modules = list(modules_needed or [])
    shell_cmds = list(shell_commands or [])

    discovered = _discover_partition_data()
    partitions = []
    for p_name in sorted(discovered.keys()):
        p = discovered[p_name]
        p.gpu_syntax = gpu_syntax
        p.use_mem_option = use_mem_option
        p.command = default_command
        p.modules_needed = list(modules)
        p.shell_commands = list(shell_cmds)
        p.description = p.description or f"Auto-discovered partition '{p.name}'"
        partitions.append(p.to_partition())

    return MachineConfig(name=cluster_name, description=desc, partitions=partitions)


def write_discovered_machine_config(
    output_path: Optional[str] = None,
    machine_name: Optional[str] = None,
    description: Optional[str] = None,
    default_command: str = DEFAULT_COMMAND,
    modules_needed: Optional[list[str]] = None,
    shell_commands: Optional[list[str]] = None,
    gpu_syntax: str = "gres",
    use_mem_option: bool = True,
) -> Path:
    """Discover machine config and write JSON to disk.

    If output_path is not provided, writes to ``~/.fleuriste/machines/<name>.json``.
    """
    config = discover_machine_config(
        machine_name=machine_name,
        description=description,
        default_command=default_command,
        modules_needed=modules_needed,
        shell_commands=shell_commands,
        gpu_syntax=gpu_syntax,
        use_mem_option=use_mem_option,
    )

    if output_path:
        out = Path(output_path).expanduser()
    else:
        out = get_machines_directory() / f"{config.name}.json"

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(config.to_dict(), handle, indent=2)
        handle.write("\n")

    return out


def _parse_csv_list(value: Optional[str]) -> list[str]:
    """Parse comma-separated list from CLI argument."""
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _build_arg_parser():
    """Create argparse parser lazily to avoid import overhead."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Discover SLURM machine properties using sinfo/scontrol and "
            "write a MachineConfig-compatible JSON file."
        )
    )
    parser.add_argument("--output", help="Output path for JSON file")
    parser.add_argument("--name", help="Machine name override")
    parser.add_argument("--description", help="Machine description override")
    parser.add_argument(
        "--command",
        default=DEFAULT_COMMAND,
        help=f"Default command written per partition (default: {DEFAULT_COMMAND})",
    )
    parser.add_argument(
        "--modules",
        help="Comma-separated modules to store in each partition",
    )
    parser.add_argument(
        "--shell-commands",
        help="Comma-separated shell commands to store in each partition",
    )
    parser.add_argument(
        "--gpu-syntax",
        choices=["gpus", "gres"],
        default="gres",
        help="GPU option syntax to store in partition configs",
    )
    parser.add_argument(
        "--no-mem-option",
        action="store_true",
        help="Store use_mem_option=false for all partitions",
    )
    return parser


def main() -> int:
    """CLI entrypoint."""
    parser = _build_arg_parser()
    args = parser.parse_args()

    try:
        out = write_discovered_machine_config(
            output_path=args.output,
            machine_name=args.name,
            description=args.description,
            default_command=args.command,
            modules_needed=_parse_csv_list(args.modules),
            shell_commands=_parse_csv_list(args.shell_commands),
            gpu_syntax=args.gpu_syntax,
            use_mem_option=not args.no_mem_option,
        )
    except SlurmDiscoveryError as exc:
        print(f"Discovery failed: {exc}")
        return 2

    print(f"Wrote machine config: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

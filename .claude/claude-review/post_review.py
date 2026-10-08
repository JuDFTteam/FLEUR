#!/usr/bin/env python3
"""
Post a Claude Code review to a GitLab merge request.

  - findings that land on a line visible in the MR diff -> inline diff discussions
  - overall summary + findings that can't be anchored   -> one summary note (updated in place)

Environment (GitLab sets the CI_* ones in merge request pipelines):
  CI_API_V4_URL, CI_PROJECT_ID, CI_MERGE_REQUEST_IID
  GITLAB_REVIEW_TOKEN          project access token, scope "api", role Reporter+
  BASE_SHA, START_SHA, HEAD_SHA  from the MR "versions" API (see .gitlab-ci.yml)

Usage: post_review.py claude_out.json
No third-party dependencies (stdlib only).
"""
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request

API = os.environ["CI_API_V4_URL"]
PROJECT = os.environ["CI_PROJECT_ID"]
MR_IID = os.environ["CI_MERGE_REQUEST_IID"]
TOKEN = os.environ["GITLAB_REVIEW_TOKEN"]
BASE, START, HEAD = (os.environ[k] for k in ("BASE_SHA", "START_SHA", "HEAD_SHA"))

MARKER = "<!-- claude-review -->"
ORDER = ["blocker", "major", "minor", "nit"]
ICON = {"blocker": "🔴", "major": "🟠", "minor": "🟡", "nit": "⚪"}


# --------------------------------------------------------------------------- GitLab API
def api(method, path, data=None, params=None):
    url = f"{API}/projects/{PROJECT}/merge_requests/{MR_IID}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url,
        method=method,
        data=None if data is None else json.dumps(data).encode(),
        headers={"PRIVATE-TOKEN": TOKEN, "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp), resp.headers


def api_all(path):
    items, page = [], 1
    while page:
        data, headers = api("GET", path, params={"per_page": 100, "page": page})
        items.extend(data)
        page = int(headers.get("X-Next-Page") or 0)
    return items


# --------------------------------------------------------------------------- diff parsing
def diff_line_map():
    """
    Lines GitLab can attach a diff comment to.
    Returns {new_path: {"old_path": str, "lines": {new_line: old_line | None}}}
    old_line is None for added lines, set for unchanged context lines
    (GitLab needs both numbers for context lines).
    """
    diff = subprocess.run(
        ["git", "diff", "--no-color", "--find-renames", "-U3", BASE, HEAD],
        capture_output=True, text=True, check=True,
    ).stdout

    files, cur, old_path, in_hunk = {}, None, None, False
    o = n = 0
    for line in diff.splitlines():
        if line.startswith("diff --git "):
            cur, old_path, in_hunk = None, None, False
        elif not in_hunk and line.startswith("--- "):
            old_path = None if line == "--- /dev/null" else line[6:]
        elif not in_hunk and line.startswith("+++ "):
            if line != "+++ /dev/null":  # deleted file -> nothing to comment on
                new_path = line[6:]
                cur = files.setdefault(new_path, {"old_path": old_path or new_path, "lines": {}})
        elif line.startswith("@@"):
            m = re.match(r"@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@", line)
            o, n, in_hunk = int(m[1]), int(m[2]), True
        elif in_hunk and cur is not None:
            tag = line[:1]
            if tag == " ":
                cur["lines"][n] = o
                o += 1
                n += 1
            elif tag == "+":
                cur["lines"][n] = None
                n += 1
            elif tag == "-":
                o += 1
    return files


# --------------------------------------------------------------------------- Claude output
def load_review(path):
    """Extract {"summary", "findings"} from `claude -p --output-format json` output."""
    try:
        with open(path) as fh:
            out = json.load(fh)
    except (OSError, ValueError) as e:
        return {"summary": f"⚠️ Claude run produced no usable output ({e}).", "findings": []}, None

    text = out.get("result") or ""
    cost = out.get("total_cost_usd")
    match = re.search(r"\{.*\}", text, re.S)  # tolerate stray prose / code fences
    try:
        review = json.loads(match.group(0))
        review.setdefault("findings", [])
    except (AttributeError, ValueError):
        raw = text or f"(no result; run ended with subtype `{out.get('subtype')}`)"
        review = {"summary": "⚠️ Could not parse a structured review. Raw output:\n\n" + raw,
                  "findings": []}
    return review, cost


def severity_rank(f):
    s = f.get("severity")
    return ORDER.index(s) if s in ORDER else len(ORDER)


def finding_body(f):
    sev = f.get("severity", "minor")
    body = f"{MARKER}\n{ICON.get(sev, '🔹')} **{sev}**: {str(f.get('comment', '')).strip()}"
    if f.get("suggestion"):
        body += f"\n\n```suggestion\n{f['suggestion']}\n```"
    return body


# --------------------------------------------------------------------------- main
def main():
    review, cost = load_review(sys.argv[1])
    findings = sorted(review["findings"], key=severity_rank)
    diff = diff_line_map()

    # Inline comments posted by earlier runs (so a new push doesn't duplicate them)
    existing = set()
    for d in api_all("/discussions"):
        note = d["notes"][0]
        pos = note.get("position")
        if pos and MARKER in note["body"]:
            existing.add((pos.get("new_path"), pos.get("new_line")))

    unanchored, posted, skipped = [], 0, 0
    for f in findings:
        path, line = f.get("file"), f.get("line")
        entry = diff.get(path)
        if not isinstance(line, int) or not entry or line not in entry["lines"]:
            unanchored.append(f)
            continue
        if (path, line) in existing:
            skipped += 1
            continue

        position = {
            "position_type": "text",
            "base_sha": BASE, "start_sha": START, "head_sha": HEAD,
            "old_path": entry["old_path"], "new_path": path,
            "new_line": line,
        }
        if entry["lines"][line] is not None:  # unchanged context line
            position["old_line"] = entry["lines"][line]

        try:
            api("POST", "/discussions", {"body": finding_body(f), "position": position})
            posted += 1
        except urllib.error.HTTPError as e:
            msg = e.read().decode(errors="replace")[:300]
            print(f"inline comment on {path}:{line} rejected ({e.code}: {msg}); "
                  "moving it to the summary", file=sys.stderr)
            unanchored.append(f)

    # ----- summary note
    out = [MARKER, "## 🤖 Claude review", "", str(review.get("summary", "")).strip(), ""]
    if findings:
        counts = {}
        for f in findings:
            s = f.get("severity", "minor")
            counts[s] = counts.get(s, 0) + 1
        tally = ", ".join(f"{ICON.get(s, '🔹')} {c} {s}" for s, c in counts.items())
        out.append(f"**Findings:** {tally} · {posted} posted inline"
                   + (f" · {skipped} already posted earlier" if skipped else ""))
    else:
        out.append("No issues found.")

    if unanchored:
        out += ["", "### Findings outside the diff", ""]
        for f in unanchored:
            sev = f.get("severity", "minor")
            loc = f"`{f['file']}:{f.get('line', '?')}`" if f.get("file") else "general"
            comment = str(f.get("comment", "")).strip().replace("\n", "\n  ")
            out.append(f"- {ICON.get(sev, '🔹')} **{sev}** {loc}: {comment}")

    footer = f"Reviewed commit `{HEAD[:8]}`"
    if cost:
        footer += f" · ≈ ${cost:.2f} (API-equivalent)"
    out += ["", f"<sub>{footer}. Automated review, please verify before acting.</sub>"]
    body = "\n".join(out)

    old = next((nt for nt in api_all("/notes")
                if MARKER in nt["body"] and nt.get("type") is None), None)
    if old:
        api("PUT", f"/notes/{old['id']}", {"body": body})
    else:
        api("POST", "/notes", {"body": body})

    print(f"posted {posted} inline, {len(unanchored)} in summary, {skipped} skipped as duplicates")


if __name__ == "__main__":
    main()

# FLEURiste MCP Server

A containerised [Model Context Protocol (MCP)](https://modelcontextprotocol.io/)
server that exposes FLEURiste's FLEUR inp.xml editing capabilities as tools for
LLM-based agents.

## Quick start

```bash
cd packaging/docker/fleuriste-mcp

# Build and launch
podman-compose up --build -d

# Check it's running
curl -s http://localhost:8000/mcp | head

# Follow logs
podman-compose logs -f

# Stop
podman-compose down
```

## MCP endpoint

| Transport | URL |
|---|---|
| Streamable HTTP (default) | `http://localhost:8000/mcp` |
| SSE (set `MCP_TRANSPORT=sse`) | `http://localhost:8000/sse` |

## Working with FLEUR files

Mount your calculation directory into `/work`:

```bash
# Using compose override
podman-compose run -v /path/to/calculation:/work fleuriste-mcp

# Or directly
podman run --rm -p 8000:8000 \
  -v $PWD:/work:z \
  fleuriste-mcp:latest
```

## Available MCP tools

| Tool | Description |
|---|---|
| `schema_search` | Search the FLEUR XML schema for elements, attributes, documentation |
| `schema_element_info` | Detailed info for a schema element path |
| `read_inp_xml` | Read an inp.xml file |
| `get_xml_tree` | Compact tree overview of the XML structure |
| `get_element` | Get XML subtree at a given path |
| `set_attribute` | Set an attribute value on element(s) |
| `set_text` | Set text content of element(s) |
| `add_element` | Add a new child element (schema-aware) |
| `remove_element` | Remove an element |
| `validate` | Validate inp.xml against the FLEUR schema |
| `save_inp_xml` | Save modifications back to disk |
| `kpoints_list` | List k-point sets in inp.xml |
| `kpoints_show` | Show details of a k-point list |
| `kpoints_select` | Set the active k-point list |
| `kpoints_delete` | Delete a k-point list |
| `analyze_fleur_input` | Analyse calculation setup (atoms, memory, matrix size…) |
| `list_workdir` | List files in /work |
| `read_file_content` | Read any text file from /work |

## Configuration

Environment variables (set in `compose.yaml` or on the command line):

| Variable | Default | Description |
|---|---|---|
| `MCP_TRANSPORT` | `streamable-http` | `streamable-http` or `sse` |
| `MCP_HOST` | `0.0.0.0` | Bind address |
| `MCP_PORT` | `8000` | Bind port |
| `FLEURISTE_SCHEMA` | `/app/FleurInputSchema.xsd` | Path to the XSD schema |
| `FLEURISTE_WORKDIR` | `/work` | Base directory for file operations |

## Connecting from an MCP client

### VS Code / Copilot

Add to your `.vscode/mcp.json`:

```json
{
  "servers": {
    "fleuriste": {
      "type": "http",
      "url": "http://localhost:8000/mcp"
    }
  }
}
```

### Claude Desktop

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "fleuriste": {
      "url": "http://localhost:8000/mcp"
    }
  }
}
```

### Python (mcp client library)

```python
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client

async with streamablehttp_client("http://localhost:8000/mcp") as (r, w, _):
    async with ClientSession(r, w) as session:
        await session.initialize()
        tools = await session.list_tools()
        print([t.name for t in tools.tools])

        result = await session.call_tool("schema_search", {"query": "cutoffs"})
        print(result.content[0].text)
```

## Building without compose

```bash
# From the repo root:
podman build -t fleuriste-mcp:latest \
  -f packaging/docker/fleuriste-mcp/Containerfile .
```

## Files

```
packaging/docker/fleuriste-mcp/
├── Containerfile      # Multi-stage container build
├── compose.yaml       # podman-compose service definition
├── mcp_server.py      # MCP server (FastMCP + fleuriste library)
└── README.md          # This file
```

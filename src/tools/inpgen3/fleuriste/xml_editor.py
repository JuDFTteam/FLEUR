"""
XML Document Model for FLEUR input files.

Provides a tree-based model for XML documents with schema-aware
editing capabilities.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Iterator, Tuple
from pathlib import Path
import xml.etree.ElementTree as ET
from xml.dom import minidom
import copy
import re
import itertools

from .schema_parser import XSDParser, SchemaElement, AttributeInfo


# Common namespace mappings
KNOWN_NAMESPACES = {
    "http://www.w3.org/2001/XInclude": "xi",
    "http://www.w3.org/2001/XMLSchema-instance": "xsi",
}

# Global counter for unique node IDs (avoids id() reuse issues with large trees)
_node_id_counter = itertools.count(1)


@dataclass
class XMLNode:
    """Represents a node in the XML document tree."""
    tag: str
    attributes: Dict[str, str] = field(default_factory=dict)
    text: Optional[str] = None
    tail: Optional[str] = None
    children: List["XMLNode"] = field(default_factory=list)
    parent: Optional["XMLNode"] = None
    schema_element: Optional[SchemaElement] = None
    # Store namespace prefix mappings for this element
    namespaces: Dict[str, str] = field(default_factory=dict)
    
    # Unique identifier for this node instance (uses global counter to avoid id reuse)
    _id: int = field(default_factory=lambda: next(_node_id_counter))
    
    @property
    def uid(self) -> int:
        """Unique identifier for this node."""
        return self._id
    
    def add_child(self, child: "XMLNode") -> "XMLNode":
        """Add a child node."""
        child.parent = self
        self.children.append(child)
        return child
    
    def remove_child(self, child: "XMLNode"):
        """Remove a child node."""
        if child in self.children:
            self.children.remove(child)
            child.parent = None
    
    def insert_child(self, index: int, child: "XMLNode"):
        """Insert a child at a specific position."""
        child.parent = self
        self.children.insert(index, child)
    
    def get_child_by_tag(self, tag: str) -> Optional["XMLNode"]:
        """Get first child with given tag."""
        for child in self.children:
            if child.tag == tag:
                return child
        return None
    
    def get_children_by_tag(self, tag: str) -> List["XMLNode"]:
        """Get all children with given tag."""
        return [c for c in self.children if c.tag == tag]
    
    def _get_readable_tag(self) -> str:
        """Get readable tag name, converting Clark notation to prefix:localname."""
        if self.tag.startswith("{"):
            ns_end = self.tag.find("}")
            if ns_end > 0:
                ns_uri = self.tag[1:ns_end]
                local_name = self.tag[ns_end + 1:]
                prefix = self.namespaces.get(ns_uri) or KNOWN_NAMESPACES.get(ns_uri)
                if prefix:
                    return f"{prefix}:{local_name}"
                return local_name
        return self.tag
    
    @property
    def display_tag(self) -> str:
        """Get human-readable tag name (resolves namespace prefixes)."""
        return self._get_readable_tag()
    
    def get_path(self) -> List[str]:
        """Get the path from root to this node."""
        path = []
        node = self
        while node is not None:
            path.insert(0, node._get_readable_tag())
            node = node.parent
        return path
    
    def get_path_string(self) -> str:
        """Get path as a string."""
        return "/" + "/".join(self.get_path())
    
    def clone(self) -> "XMLNode":
        """Create a deep copy of this node."""
        new_node = XMLNode(
            tag=self.tag,
            attributes=dict(self.attributes),
            text=self.text,
            tail=self.tail,
            namespaces=dict(self.namespaces)
        )
        for child in self.children:
            cloned_child = child.clone()
            new_node.add_child(cloned_child)
        return new_node
    
    def _get_prefixed_tag(self, tag: str, ns_map: Dict[str, str]) -> str:
        """Convert {uri}localname tag to prefix:localname format."""
        if tag.startswith("{"):
            ns_end = tag.find("}")
            if ns_end > 0:
                ns_uri = tag[1:ns_end]
                local_name = tag[ns_end + 1:]
                # Look up prefix in ns_map or known namespaces
                prefix = ns_map.get(ns_uri) or KNOWN_NAMESPACES.get(ns_uri)
                if prefix:
                    return f"{prefix}:{local_name}"
        return tag
    
    def _collect_namespaces(self) -> Dict[str, str]:
        """Collect all namespace URIs used in this node and descendants."""
        ns_map = dict(self.namespaces)
        for child in self.children:
            ns_map.update(child._collect_namespaces())
        return ns_map
    
    def to_element(self, ns_map: Optional[Dict[str, str]] = None) -> ET.Element:
        """Convert to an ElementTree Element."""
        if ns_map is None:
            ns_map = {}
        
        # Merge with our local namespaces
        combined_ns = dict(ns_map)
        combined_ns.update(self.namespaces)
        
        # Get the prefixed tag name
        prefixed_tag = self._get_prefixed_tag(self.tag, combined_ns)
        
        # Build attributes, handling namespace prefixes
        new_attribs = {}
        for key, value in self.attributes.items():
            prefixed_key = self._get_prefixed_tag(key, combined_ns)
            new_attribs[prefixed_key] = value
        
        # Add namespace declarations
        for ns_uri, prefix in self.namespaces.items():
            new_attribs[f"xmlns:{prefix}"] = ns_uri
        
        elem = ET.Element(prefixed_tag, new_attribs)
        elem.text = self.text
        elem.tail = self.tail
        for child in self.children:
            elem.append(child.to_element(combined_ns))
        return elem
    
    @classmethod
    def from_element(cls, elem: ET.Element, parent: Optional["XMLNode"] = None, 
                     inherited_ns: Optional[Dict[str, str]] = None) -> "XMLNode":
        """Create XMLNode from an ElementTree Element."""
        if inherited_ns is None:
            inherited_ns = {}
        
        # Extract namespace declarations from attributes (rarely present with ET)
        namespaces = {}
        regular_attribs = {}
        
        for key, value in elem.attrib.items():
            if key.startswith("{"):
                # Already a Clark notation attribute (e.g., from namespace)
                regular_attribs[key] = value
            elif key.startswith("xmlns:"):
                # Namespace declaration: xmlns:prefix="uri"
                prefix = key[6:]  # Remove "xmlns:"
                namespaces[value] = prefix
            elif key == "xmlns":
                # Default namespace declaration
                namespaces[value] = ""
            else:
                regular_attribs[key] = value
        
        # Detect namespace from tag in Clark notation {uri}localname
        tag = elem.tag
        if tag.startswith("{"):
            ns_end = tag.find("}")
            if ns_end > 0:
                ns_uri = tag[1:ns_end]
                # Check if this namespace is already known (inherited or in KNOWN_NAMESPACES)
                if ns_uri not in inherited_ns and ns_uri not in namespaces:
                    # Try to find a known prefix for this namespace
                    prefix = KNOWN_NAMESPACES.get(ns_uri)
                    if prefix:
                        namespaces[ns_uri] = prefix
        
        # Also check attributes for namespaced attributes
        for key in list(regular_attribs.keys()):
            if key.startswith("{"):
                ns_end = key.find("}")
                if ns_end > 0:
                    ns_uri = key[1:ns_end]
                    if ns_uri not in inherited_ns and ns_uri not in namespaces:
                        prefix = KNOWN_NAMESPACES.get(ns_uri)
                        if prefix:
                            namespaces[ns_uri] = prefix
        
        # Merge inherited namespaces for children
        combined_ns = dict(inherited_ns)
        combined_ns.update(namespaces)
        
        node = cls(
            tag=elem.tag,
            attributes=regular_attribs,
            text=elem.text.strip() if elem.text and elem.text.strip() else None,
            tail=elem.tail.strip() if elem.tail and elem.tail.strip() else None,
            parent=parent,
            namespaces=namespaces
        )
        for child_elem in elem:
            child_node = cls.from_element(child_elem, parent=node, inherited_ns=combined_ns)
            node.children.append(child_node)
        return node
    
    def iter_descendants(self) -> Iterator["XMLNode"]:
        """Iterate over all descendant nodes."""
        for child in self.children:
            yield child
            yield from child.iter_descendants()
    
    def find_by_uid(self, uid: int) -> Optional["XMLNode"]:
        """Find a node by its unique ID."""
        if self._id == uid:
            return self
        for child in self.children:
            result = child.find_by_uid(uid)
            if result:
                return result
        return None
    
    def get_display_name(self) -> str:
        """Get a display name for this node."""
        display_tag = self._get_readable_tag()
        
        # Try to use meaningful attributes for identification
        name_attrs = ["name", "species", "element", "label", "listName", "href"]
        for attr in name_attrs:
            if attr in self.attributes:
                return f"{display_tag} ({self.attributes[attr]})"
        return display_tag
    
    def has_content(self) -> bool:
        """Check if node has text content."""
        return self.text is not None and len(self.text.strip()) > 0


def _process_xincludes(elem: ET.Element, base_dir: Path) -> None:
    """Process XInclude directives in-place.

    Each ``<xi:include href="..."/>`` element whose referenced file exists is
    replaced with the parsed root of that file.  If the file cannot be found
    (or the href is a remote URL / text include) the xi:include element is left
    untouched so that round-trip serialisation preserves it.
    """
    XI_INCLUDE = "{http://www.w3.org/2001/XInclude}include"

    i = 0
    while i < len(elem):
        child = elem[i]
        if child.tag == XI_INCLUDE:
            href = child.get("href", "")
            parse = child.get("parse", "xml")
            replaced = False
            # Only handle local XML includes (skip URLs and text-mode includes)
            if href and parse == "xml" and not href.startswith(("http://", "https://", "ftp://")):
                include_path = base_dir / href
                if include_path.exists():
                    try:
                        included_root = ET.parse(str(include_path)).getroot()
                        # Preserve the whitespace tail that was on the xi:include tag
                        included_root.tail = child.tail
                        elem.remove(child)
                        elem.insert(i, included_root)
                        # Recurse so nested xi:includes are also resolved
                        _process_xincludes(included_root, include_path.parent)
                        replaced = True
                    except ET.ParseError:
                        pass  # keep the xi:include element on parse failure
            if not replaced:
                i += 1
        else:
            _process_xincludes(child, base_dir)
            i += 1


class XMLDocument:
    """Represents a complete XML document with schema awareness."""
    
    def __init__(self, schema: Optional[XSDParser] = None):
        self.schema = schema
        self.root: Optional[XMLNode] = None
        self.file_path: Optional[Path] = None
        self._modified = False
        self._undo_stack: List[XMLNode] = []
        self._redo_stack: List[XMLNode] = []
    
    @property
    def is_modified(self) -> bool:
        return self._modified
    
    def mark_modified(self):
        """Mark document as modified."""
        self._modified = True
    
    def clear_modified(self):
        """Clear modified flag."""
        self._modified = False
    
    def load(self, file_path: str | Path):
        """Load an XML document from file."""
        self.file_path = Path(file_path)
        tree = ET.parse(self.file_path)
        # Resolve xi:include directives; missing files are left as xi:include nodes
        _process_xincludes(tree.getroot(), self.file_path.parent)
        self.root = XMLNode.from_element(tree.getroot())
        self._link_schema()
        self._modified = False
    
    def load_string(self, xml_string: str):
        """Load an XML document from a string."""
        root_elem = ET.fromstring(xml_string)
        self.root = XMLNode.from_element(root_elem)
        self._link_schema()
        self._modified = False
    
    def save(self, file_path: Optional[str | Path] = None):
        """Save the document to file."""
        path = Path(file_path) if file_path else self.file_path
        if not path:
            raise ValueError("No file path specified")
        
        if not self.root:
            raise ValueError("No document to save")
        
        # Convert to ElementTree and save
        root_elem = self.root.to_element()
        
        # Pretty print
        xml_string = ET.tostring(root_elem, encoding="unicode")
        dom = minidom.parseString(xml_string)
        pretty_xml = dom.toprettyxml(indent="   ")
        
        # Remove extra blank lines and fix declaration
        lines = pretty_xml.split("\n")
        lines = [l for l in lines if l.strip()]
        # Update XML declaration
        if lines and lines[0].startswith("<?xml"):
            lines[0] = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
        
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        
        self.file_path = path
        self._modified = False
    
    def to_string(self) -> str:
        """Convert document to XML string."""
        if not self.root:
            return ""
        
        root_elem = self.root.to_element()
        xml_string = ET.tostring(root_elem, encoding="unicode")
        dom = minidom.parseString(xml_string)
        return dom.toprettyxml(indent="   ")
    
    def _link_schema(self):
        """Link schema elements to document nodes."""
        if not self.schema or not self.root:
            return
        
        root_schema = self.schema.get_root_element()
        if root_schema:
            self._link_schema_recursive(self.root, root_schema)
    
    def _link_schema_recursive(self, node: XMLNode, schema_elem: SchemaElement):
        """Recursively link schema elements to nodes."""
        node.schema_element = schema_elem
        
        for child_node in node.children:
            # Get the local name (without namespace) for matching
            child_tag = child_node.tag
            if child_tag.startswith("{"):
                ns_end = child_tag.find("}")
                if ns_end > 0:
                    child_tag = child_tag[ns_end + 1:]
            
            # Find matching schema element
            for child_schema in schema_elem.children:
                if child_schema.name == child_tag or child_schema.name == child_node.tag:
                    self._link_schema_recursive(child_node, child_schema)
                    break
    
    def save_undo_state(self):
        """Save current state for undo."""
        if self.root:
            self._undo_stack.append(self.root.clone())
            self._redo_stack.clear()
            # Limit undo stack size
            if len(self._undo_stack) > 50:
                self._undo_stack.pop(0)
    
    def undo(self) -> bool:
        """Undo last change."""
        if self._undo_stack and self.root:
            self._redo_stack.append(self.root.clone())
            self.root = self._undo_stack.pop()
            self._link_schema()
            return True
        return False
    
    def redo(self) -> bool:
        """Redo last undone change."""
        if self._redo_stack and self.root:
            self._undo_stack.append(self.root.clone())
            self.root = self._redo_stack.pop()
            self._link_schema()
            return True
        return False
    
    def get_addable_children(self, node: XMLNode) -> List[SchemaElement]:
        """Get list of elements that can be added as children of the given node."""
        if not node.schema_element:
            return []
        
        addable = []
        for child_schema in node.schema_element.children:
            # Check if element can be added (respecting maxOccurs)
            existing_count = len(node.get_children_by_tag(child_schema.name))
            if child_schema.max_occurs == -1 or existing_count < child_schema.max_occurs:
                addable.append(child_schema)
        
        return addable
    
    def get_required_children(self, node: XMLNode) -> List[SchemaElement]:
        """Get list of required children that are missing."""
        if not node.schema_element:
            return []
        
        missing = []
        for child_schema in node.schema_element.children:
            if child_schema.is_required:
                if not node.get_child_by_tag(child_schema.name):
                    missing.append(child_schema)
        
        return missing
    
    def create_element_from_schema(self, schema_elem: SchemaElement) -> XMLNode:
        """Create a new XML node from a schema element with defaults."""
        node = XMLNode(tag=schema_elem.name)
        
        # Add required attributes with defaults
        for attr_name, attr_info in schema_elem.attributes.items():
            if attr_info.is_required:
                default_val = attr_info.default or self._get_default_for_type(attr_info.type_name)
                node.attributes[attr_name] = default_val
            elif attr_info.default:
                node.attributes[attr_name] = attr_info.default
        
        # Add simple content if needed
        if schema_elem.is_simple_content:
            node.text = self._get_default_for_type(schema_elem.simple_content_type or "")
        
        # Recursively add required children
        for child_schema in schema_elem.children:
            if child_schema.is_required:
                child_node = self.create_element_from_schema(child_schema)
                node.add_child(child_node)
        
        node.schema_element = schema_elem
        return node
    
    def _get_default_for_type(self, type_name: str) -> str:
        """Get a sensible default value for a type."""
        if not type_name:
            return ""
        
        if self.schema and self.schema.is_boolean_type(type_name):
            return "F"
        
        if self.schema and self.schema.is_numeric_type(type_name):
            return "0"
        
        if self.schema:
            enums = self.schema.get_enum_values(type_name)
            if enums:
                return enums[0]
        
        return ""
    
    def validate_node(self, node: XMLNode) -> List[str]:
        """Validate a node against the schema. Returns list of error messages."""
        errors = []
        
        if not node.schema_element:
            return errors
        
        schema = node.schema_element
        
        # Check required attributes
        for attr_name, attr_info in schema.attributes.items():
            if attr_info.is_required and attr_name not in node.attributes:
                errors.append(f"Missing required attribute: {attr_name}")
        
        # Check required children
        for child_schema in schema.children:
            if child_schema.is_required:
                if not node.get_child_by_tag(child_schema.name):
                    errors.append(f"Missing required element: {child_schema.name}")
        
        # Check for unexpected attributes
        if schema.attributes:
            for attr_name in node.attributes:
                if attr_name not in schema.attributes:
                    errors.append(f"Unexpected attribute: {attr_name}")
        
        return errors
    
    def find_nodes_by_tag(self, tag: str) -> List[XMLNode]:
        """Find all nodes with the given tag."""
        if not self.root:
            return []
        
        results = []
        if self.root.tag == tag:
            results.append(self.root)
        
        for node in self.root.iter_descendants():
            if node.tag == tag:
                results.append(node)
        
        return results

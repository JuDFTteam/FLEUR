"""
XSD Schema Parser for FLEUR input files.

Parses the FleurInputSchema.xsd and provides structured access to
element definitions, types, and constraints.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Set
from pathlib import Path
import xml.etree.ElementTree as ET


XSD_NS = "http://www.w3.org/2001/XMLSchema"
XSD = f"{{{XSD_NS}}}"


@dataclass
class SearchResult:
    """A search result with path and match information."""
    path: str  # XML path like /fleurInput/calculationSetup/cutoffs
    element_name: str
    match_type: str  # "element", "attribute", "documentation"
    match_text: str  # The text that matched
    context: str  # Additional context (e.g., documentation snippet)


@dataclass
class AttributeInfo:
    """Information about an XML attribute from the schema."""
    name: str
    type_name: str
    use: str = "optional"  # "required" or "optional"
    default: Optional[str] = None
    documentation: str = ""
    
    @property
    def is_required(self) -> bool:
        return self.use == "required"


@dataclass
class SchemaElement:
    """Represents an element definition from the XSD schema."""
    name: str
    type_name: Optional[str] = None
    min_occurs: int = 1
    max_occurs: int = 1  # -1 means unbounded
    attributes: Dict[str, AttributeInfo] = field(default_factory=dict)
    children: List["SchemaElement"] = field(default_factory=list)
    choice_groups: List[List[str]] = field(default_factory=list)  # Groups of mutually exclusive elements
    is_simple_content: bool = False
    simple_content_type: Optional[str] = None
    documentation: str = ""
    
    @property
    def is_optional(self) -> bool:
        return self.min_occurs == 0
    
    @property
    def is_unbounded(self) -> bool:
        return self.max_occurs == -1
    
    @property
    def is_required(self) -> bool:
        return self.min_occurs > 0


@dataclass
class EnumType:
    """Represents an enumeration type from the schema."""
    name: str
    values: List[str]


@dataclass
class SimpleType:
    """Represents a simple type definition."""
    name: str
    base_type: str
    pattern: Optional[str] = None
    enumerations: List[str] = field(default_factory=list)
    min_value: Optional[str] = None
    max_value: Optional[str] = None
    length: Optional[int] = None


class XSDParser:
    """Parser for FLEUR XSD schema files."""
    
    def __init__(self, schema_path: str | Path):
        self.schema_path = Path(schema_path)
        self.tree = ET.parse(self.schema_path)
        self.root = self.tree.getroot()
        
        # Storage for parsed types
        self.complex_types: Dict[str, ET.Element] = {}
        self.simple_types: Dict[str, SimpleType] = {}
        self.enum_types: Dict[str, EnumType] = {}
        self.elements: Dict[str, SchemaElement] = {}
        
        # Parse the schema
        self._parse_schema()
    
    def _parse_schema(self):
        """Parse all type definitions from the schema."""
        # First pass: collect all type definitions (both simple and complex)
        root_element_def = None
        for child in self.root:
            tag = child.tag.replace(XSD, "")
            name = child.get("name", "")
            
            if tag == "complexType":
                self.complex_types[name] = child
            elif tag == "simpleType":
                self._parse_simple_type(child, name)
            elif tag == "element":
                # Store root element for second pass
                root_element_def = child
        
        # Second pass: parse root element (now all types are available)
        if root_element_def is not None:
            self._parse_root_element(root_element_def)
    
    def _parse_simple_type(self, elem: ET.Element, name: str):
        """Parse a simpleType definition."""
        restriction = elem.find(f"{XSD}restriction")
        union = elem.find(f"{XSD}union")
        list_elem = elem.find(f"{XSD}list")
        
        if restriction is not None:
            base_type = restriction.get("base", "xsd:string")
            simple_type = SimpleType(name=name, base_type=base_type)
            
            # Check for pattern
            pattern = restriction.find(f"{XSD}pattern")
            if pattern is not None:
                simple_type.pattern = pattern.get("value")
            
            # Check for enumerations
            enums = restriction.findall(f"{XSD}enumeration")
            if enums:
                simple_type.enumerations = [e.get("value", "") for e in enums]
                self.enum_types[name] = EnumType(name=name, values=simple_type.enumerations)
            
            # Check for min/max
            min_incl = restriction.find(f"{XSD}minInclusive")
            max_incl = restriction.find(f"{XSD}maxInclusive")
            if min_incl is not None:
                simple_type.min_value = min_incl.get("value")
            if max_incl is not None:
                simple_type.max_value = max_incl.get("value")
            
            # Check for length
            length = restriction.find(f"{XSD}length")
            if length is not None:
                simple_type.length = int(length.get("value", "0"))
            
            self.simple_types[name] = simple_type
            
        elif union is not None:
            # Union type - combine member types
            member_types = union.get("memberTypes", "").split()
            simple_type = SimpleType(name=name, base_type="union")
            self.simple_types[name] = simple_type
            
        elif list_elem is not None:
            # List type
            item_type = list_elem.get("itemType", "xsd:string")
            simple_type = SimpleType(name=name, base_type=f"list({item_type})")
            self.simple_types[name] = simple_type
    
    def _parse_root_element(self, elem: ET.Element):
        """Parse the root element definition."""
        name = elem.get("name", "")
        type_name = elem.get("type")
        
        # Extract documentation
        documentation = self._get_documentation(elem)
        
        root_elem = SchemaElement(name=name, type_name=type_name, documentation=documentation)
        
        if type_name and type_name in self.complex_types:
            type_elem = self.complex_types[type_name]
            self._populate_element_from_type(root_elem, type_elem)
            # Also get documentation from the type if element doesn't have its own
            if not root_elem.documentation:
                root_elem.documentation = self._get_documentation(type_elem)
        
        self.elements[name] = root_elem
    
    def _populate_element_from_type(self, schema_elem: SchemaElement, type_elem: ET.Element):
        """Populate a SchemaElement with info from a complexType definition."""
        # Parse attributes
        for attr in type_elem.findall(f"{XSD}attribute"):
            attr_info = self._parse_attribute(attr)
            schema_elem.attributes[attr_info.name] = attr_info
        
        # Check for simpleContent
        simple_content = type_elem.find(f"{XSD}simpleContent")
        if simple_content is not None:
            extension = simple_content.find(f"{XSD}extension")
            if extension is not None:
                schema_elem.is_simple_content = True
                schema_elem.simple_content_type = extension.get("base")
                # Parse attributes from extension
                for attr in extension.findall(f"{XSD}attribute"):
                    attr_info = self._parse_attribute(attr)
                    schema_elem.attributes[attr_info.name] = attr_info
            return
        
        # Parse child elements from sequence, all, or choice
        sequence = type_elem.find(f"{XSD}sequence")
        all_elem = type_elem.find(f"{XSD}all")
        choice = type_elem.find(f"{XSD}choice")
        
        container = sequence or all_elem or choice
        if container is not None:
            self._parse_element_container(schema_elem, container)
    
    def _get_documentation(self, elem: ET.Element) -> str:
        """Extract documentation text from an element's annotation."""
        annotation = elem.find(f"{XSD}annotation")
        if annotation is not None:
            doc = annotation.find(f"{XSD}documentation")
            if doc is not None and doc.text:
                return doc.text.strip()
        return ""
    
    def _parse_attribute(self, attr_elem: ET.Element) -> AttributeInfo:
        """Parse an attribute definition."""
        return AttributeInfo(
            name=attr_elem.get("name", ""),
            type_name=attr_elem.get("type", "xsd:string"),
            use=attr_elem.get("use", "optional"),
            default=attr_elem.get("default"),
            documentation=self._get_documentation(attr_elem)
        )
    
    def _parse_element_container(self, parent: SchemaElement, container: ET.Element):
        """Parse elements from a sequence, all, or choice container."""
        for child in container:
            tag = child.tag.replace(XSD, "")
            
            if tag == "element":
                child_elem = self._parse_child_element(child)
                parent.children.append(child_elem)
            elif tag == "choice":
                choice_elements = []
                self._parse_choice_group(parent, child, choice_elements)
                if choice_elements:
                    parent.choice_groups.append(choice_elements)
            elif tag == "sequence":
                self._parse_element_container(parent, child)
            elif tag == "all":
                self._parse_element_container(parent, child)
    
    def _parse_choice_group(self, parent: SchemaElement, choice_elem: ET.Element, choice_names: List[str]):
        """Parse a choice group and add elements to parent."""
        for child in choice_elem:
            tag = child.tag.replace(XSD, "")
            if tag == "element":
                child_elem = self._parse_child_element(child)
                parent.children.append(child_elem)
                choice_names.append(child_elem.name)
    
    def _parse_child_element(self, elem: ET.Element) -> SchemaElement:
        """Parse a child element definition."""
        name = elem.get("name", "")
        type_name = elem.get("type")
        
        min_occurs_str = elem.get("minOccurs", "1")
        max_occurs_str = elem.get("maxOccurs", "1")
        
        min_occurs = int(min_occurs_str)
        max_occurs = -1 if max_occurs_str == "unbounded" else int(max_occurs_str)
        
        # Extract documentation
        documentation = self._get_documentation(elem)
        
        child_elem = SchemaElement(
            name=name,
            type_name=type_name,
            min_occurs=min_occurs,
            max_occurs=max_occurs,
            documentation=documentation
        )
        
        # Recursively populate from type if it exists
        if type_name and type_name in self.complex_types:
            type_elem = self.complex_types[type_name]
            self._populate_element_from_type(child_elem, type_elem)
            # Also get documentation from the type if element doesn't have its own
            if not child_elem.documentation:
                child_elem.documentation = self._get_documentation(type_elem)
        
        return child_elem
    
    def get_root_element(self) -> Optional[SchemaElement]:
        """Get the root element of the schema."""
        return self.elements.get("fleurInput")
    
    def get_element_by_name(self, name: str) -> Optional[SchemaElement]:
        """Get an element by name."""
        return self.elements.get(name)
    
    def get_element_for_path(self, path: List[str]) -> Optional[SchemaElement]:
        """Get schema element for a given XPath-like path."""
        root = self.get_root_element()
        if not root or not path:
            return root
        
        current = root
        for part in path:
            found = None
            for child in current.children:
                if child.name == part:
                    found = child
                    break
            if found:
                current = found
            else:
                return None
        return current
    
    def get_enum_values(self, type_name: str) -> List[str]:
        """Get enumeration values for a type."""
        if type_name in self.enum_types:
            return self.enum_types[type_name].values
        if type_name in self.simple_types:
            return self.simple_types[type_name].enumerations
        return []
    
    def is_boolean_type(self, type_name: str) -> bool:
        """Check if a type is a boolean type."""
        return type_name in ("FleurBool", "xsd:boolean")
    
    def is_numeric_type(self, type_name: str) -> bool:
        """Check if a type is a numeric type."""
        numeric_types = {
            "xsd:integer", "xsd:nonNegativeInteger", "xsd:positiveInteger",
            "xsd:double", "xsd:float", "FleurDouble",
            "xsd:int", "xsd:long", "xsd:short"
        }
        return type_name in numeric_types
    
    def get_type_description(self, type_name: str) -> str:
        """Get a human-readable description of a type."""
        if type_name in self.enum_types:
            values = self.enum_types[type_name].values[:5]
            more = "..." if len(self.enum_types[type_name].values) > 5 else ""
            return f"One of: {', '.join(values)}{more}"
        
        if type_name in self.simple_types:
            st = self.simple_types[type_name]
            if st.pattern:
                return f"Pattern: {st.pattern}"
            if st.min_value or st.max_value:
                return f"Range: [{st.min_value or '-∞'}, {st.max_value or '∞'}]"
            return st.base_type
        
        type_descriptions = {
            "xsd:string": "Text",
            "xsd:integer": "Integer",
            "xsd:nonNegativeInteger": "Non-negative integer",
            "xsd:positiveInteger": "Positive integer",
            "xsd:double": "Decimal number",
            "xsd:float": "Decimal number",
            "FleurDouble": "Number (Fleur format)",
            "FleurBool": "T or F",
            "xsd:boolean": "true or false",
        }
        return type_descriptions.get(type_name, type_name)
    
    def search(self, query: str, max_results: int = 50) -> List[SearchResult]:
        """
        Search the schema for elements, attributes, and documentation matching the query.
        
        Args:
            query: Search string (case-insensitive)
            max_results: Maximum number of results to return
            
        Returns:
            List of SearchResult objects with matching paths
        """
        results: List[SearchResult] = []
        query_lower = query.lower()
        
        root = self.get_root_element()
        if not root:
            return results
        
        self._search_element(root, "/", query_lower, results, max_results)
        return results
    
    def _search_element(
        self, 
        elem: SchemaElement, 
        parent_path: str, 
        query: str, 
        results: List[SearchResult],
        max_results: int
    ):
        """Recursively search an element and its children."""
        if len(results) >= max_results:
            return
        
        current_path = f"{parent_path}{elem.name}"
        
        # Search element name
        if query in elem.name.lower():
            results.append(SearchResult(
                path=current_path,
                element_name=elem.name,
                match_type="element",
                match_text=elem.name,
                context=elem.documentation[:100] if elem.documentation else ""
            ))
        
        # Search element documentation
        elif elem.documentation and query in elem.documentation.lower():
            # Find the matching snippet
            idx = elem.documentation.lower().find(query)
            start = max(0, idx - 30)
            end = min(len(elem.documentation), idx + len(query) + 30)
            snippet = elem.documentation[start:end]
            if start > 0:
                snippet = "..." + snippet
            if end < len(elem.documentation):
                snippet = snippet + "..."
            
            results.append(SearchResult(
                path=current_path,
                element_name=elem.name,
                match_type="documentation",
                match_text=snippet,
                context=""
            ))
        
        # Search attributes
        for attr_name, attr_info in elem.attributes.items():
            if len(results) >= max_results:
                return
                
            if query in attr_name.lower():
                results.append(SearchResult(
                    path=f"{current_path}/@{attr_name}",
                    element_name=elem.name,
                    match_type="attribute",
                    match_text=attr_name,
                    context=attr_info.documentation[:100] if attr_info.documentation else ""
                ))
            elif attr_info.documentation and query in attr_info.documentation.lower():
                idx = attr_info.documentation.lower().find(query)
                start = max(0, idx - 30)
                end = min(len(attr_info.documentation), idx + len(query) + 30)
                snippet = attr_info.documentation[start:end]
                if start > 0:
                    snippet = "..." + snippet
                if end < len(attr_info.documentation):
                    snippet = snippet + "..."
                    
                results.append(SearchResult(
                    path=f"{current_path}/@{attr_name}",
                    element_name=elem.name,
                    match_type="attribute_doc",
                    match_text=snippet,
                    context=""
                ))
        
        # Recursively search children
        for child in elem.children:
            if len(results) >= max_results:
                return
            self._search_element(child, f"{current_path}/", query, results, max_results)


def load_schema(schema_path: str | Path) -> XSDParser:
    """Load and parse an XSD schema file."""
    return XSDParser(schema_path)

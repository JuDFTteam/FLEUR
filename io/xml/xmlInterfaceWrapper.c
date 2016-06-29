/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

/*
 * Wrapper routines for XML IO - C side
 *                                GM'16
 */

#include <stdio.h>
#include <assert.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlschemas.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

static xmlDocPtr xmlDocument;
static xmlDocPtr schemaDoc;
static xmlSchemaParserCtxtPtr schemaParserCtxt;
static xmlSchemaPtr schema;
static xmlSchemaValidCtxtPtr schemaValidCtxt;
static xmlXPathContextPtr xPathCtxt;
static xmlXPathObjectPtr xPathObj;

int initializeXMLInterface()
{
   xmlDocument = NULL;
   schemaDoc = NULL;
   schemaParserCtxt = NULL;
   schema = NULL;
   schemaValidCtxt = NULL;

   xPathCtxt = NULL;
   xPathObj = NULL;
   xmlInitParser();
   return 0;
}

int parseXMLSchema(const char* schemaFilename)
{
   schemaDoc = xmlReadFile(schemaFilename, NULL, 0);
   printf("Parsing XML Schema file: %s\n", schemaFilename);
   if (schemaDoc == NULL) 
   {
      fprintf(stderr, "Failed to parse xml schema file %s\n", schemaFilename);
      return -1;
   }
   schemaParserCtxt = xmlSchemaNewDocParserCtxt (schemaDoc);
   if (schemaParserCtxt == NULL)
   {
      fprintf(stderr, "Failed to create schemaParserCtxt from xml schema.\n");
      return -1;
   }
   schema = xmlSchemaParse(schemaParserCtxt);
   if (schema == NULL)
   {
      fprintf(stderr, "Failed to create schema from xml schema.\n");
      return -1;
   }
   schemaValidCtxt = xmlSchemaNewValidCtxt(schema);
   if (schemaValidCtxt == NULL)
   {
      fprintf(stderr, "Failed to create schemaValidCtxt from xml schema.\n");
      return -1;
   }
   xmlSchemaSetValidOptions(schemaValidCtxt,XML_SCHEMA_VAL_VC_I_CREATE);
   printf("parseXMLSchema: schemaDoc: %p\n",schemaDoc);
   printf("parseXMLSchema: schemaParserCtxt: %p\n",schemaParserCtxt);
   printf("parseXMLSchema: schema: %p\n",schema);
   printf("parseXMLSchema: schemaValidCtxt: %p\n",schemaValidCtxt);

   return 0;
}

int parseXMLDocument(const char* docFilename)
{
   xmlDocument = xmlReadFile(docFilename, NULL, 0);
   if (xmlDocument == NULL) 
   {
      fprintf(stderr, "Failed to parse xml file %s\n", docFilename);
      return -1;
   }
   printf("parseXMLDocument: xmlDocument: %p\n",xmlDocument);
   return 0;
}

int validateXMLDocument()
{

   printf("validateXMLDocument: xmlDocument: %p\n",xmlDocument);

   if (schemaValidCtxt == NULL)
   {
      fprintf(stderr, "Error: schemaValidCtxt is null in validateDocument()");
      return -1;
   }
   if (xmlDocument == NULL) 
   {
      fprintf(stderr, "Error: xmlDocument is null in validateDocument()");
      return -1;
   }
   int value;
   value = xmlSchemaValidateDoc(schemaValidCtxt, xmlDocument);
   return value;
}

int initializeXPath()
{
   if (xmlDocument == NULL) 
   {
      fprintf(stderr, "Error: xmlDocument is null in initializeXPath()");
      return -1;
   }
   xPathCtxt = xmlXPathNewContext(xmlDocument);
   return 0;
}

int getNumberOfXMLNodes(const unsigned char* xPathExpression)
{
   if (xPathCtxt == NULL) 
   {
      fprintf(stderr, "Error: xPathCtxt is null in getNumberOfNodes(...)");
      return -1;
   }
   xPathObj = xmlXPathEvalExpression(xPathExpression, xPathCtxt);
   if (xPathObj == NULL)
   {
      fprintf(stderr, "Error: xPathObj is null in getNumberOfNodes(...)");
      return -1;
   }
   xmlNodeSetPtr nodes;
   nodes = NULL;
   nodes = xPathObj->nodesetval;
   if (xPathObj == NULL)
   {
      fprintf(stderr, "Error: nodes is null in getNumberOfNodes(...)");
      return -1;
   }
   int size;
   size = (nodes) ? nodes->nodeNr : 0;
   return size;
}

extern const unsigned char* getXMLAttributeValue(const unsigned char* xPathExpression)
{
   if (xPathCtxt == NULL) 
   {
      fprintf(stderr, "Error: xPathCtxt is null in getAttributeValue(...)");
      return NULL;
   }
   xPathObj = xmlXPathEvalExpression(xPathExpression, xPathCtxt);
   if (xPathObj == NULL)
   {
      fprintf(stderr, "Error: xPathObj is null in getAttributeValue(...)");
      return NULL;
   }
   xmlNodeSetPtr nodes = xPathObj->nodesetval;
   if (nodes == NULL)
   {
      fprintf(stderr, "Error: nodes is null in getAttributeValue(...)");
      return NULL;
   }

//   return NULL;
   return xmlNodeGetContent(nodes->nodeTab[0]);

}

int freeXMLResources()
{
   xmlCleanupParser();
   if(xPathObj) xmlXPathFreeObject(xPathObj);
   if(xPathCtxt) xmlXPathFreeContext(xPathCtxt);
   if(schemaDoc) xmlFreeDoc(schemaDoc);
   if(schemaParserCtxt) xmlSchemaFreeParserCtxt(schemaParserCtxt);
   if(schema) xmlSchemaFree(schema);
   if(schemaValidCtxt) xmlSchemaFreeValidCtxt(schemaValidCtxt);
   return 0;
}

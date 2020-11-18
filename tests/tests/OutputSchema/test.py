#!/usr/bin/env python3
# load libtest
import pathlib
test_loc = pathlib.Path(__file__).parent.absolute()
import sys
sys.path.insert(0, f"{test_loc}/../../")
from libtest import TestEnv
from lxml import etree
import os

try:
   te = TestEnv()
   te.log_info("OutputSchema Test")
   te.log_info(f"Test directory: {os.path.abspath(os.path.join(te.workdir,'./../..'))}")
   for root, dirs, files in os.walk(os.path.abspath(os.path.join(te.workdir,'./../..'))):
      for file in files:
         if file == 'out.xml':
            file_path = os.path.join(root,file)
            schema_path = os.path.join(root,'FleurOutputSchema.xsd')
            te.log_info(f"Testing {file_path}")
            try:
               xmltree = etree.parse(file_path)
            except:
               te.log_error(f'File not parsable, maybe something went wrong: {file_path}')
               te.errors += 1
               continue

            if not os.path.isfile(schema_path):
               te.log_info(f"No OutputSchema for: {file_path}")
               continue

            xmlschema_doc = etree.parse(schema_path)
            xmlschema = etree.XMLSchema(xmlschema_doc)

            if not xmlschema.validate(xmltree):
               # get more information on what does not validate
               parser_on_fly = etree.XMLParser(attribute_defaults=True, schema=xmlschema, encoding='utf-8')
               outxmlfile = etree.tostring(xmltree)
               message = 'Reason is unknown'
               try:
                  tree_x = etree.fromstring(outxmlfile, parser_on_fly)
               except etree.XMLSyntaxError as msg:
                  message = msg
               te.log_error(f'Output file does not validate against the schema: {message}')
               te.errors += 1

   sys.exit(te.errors)
except Exception as e:
   te.log_error(f"Error: {e}")
   sys.exit(1)

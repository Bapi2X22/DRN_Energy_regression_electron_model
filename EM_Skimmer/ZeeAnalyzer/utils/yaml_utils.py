import yaml
from collections import OrderedDict

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
   return dumper.represent_mapping(_mapping_tag, data.iteritems())

def dict_constructor(loader, node):
   return OrderedDict(loader.construct_pairs(node))

class indent_dumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(indent_dumper, self).increase_indent(flow, False)
    
yaml.add_representer( OrderedDict , dict_representer )
yaml.add_constructor( _mapping_tag, dict_constructor ) 

def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        _mapping_tag,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)

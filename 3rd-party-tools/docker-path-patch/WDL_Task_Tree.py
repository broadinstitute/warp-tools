import WDL
import hashlib
import json
import logging

from Level import Level
log = logging.getLogger()

""" Parse a miniwdl workflow object into a tree structure of nodes (determined by WDL_Task_Tree.node_objs), usage for this class should be as follows

        tree_parser = WDL_Task_Tree()
        # Start parsing from the root workflow in parent_wdl
        tree_parser.parse_root_workflow(parent_wdl.workflow)

    tree_parser.tree is a list of the parsed nodes with parent relationships described by parent_digest
    tree_parser.tree_map and tree_parser.caller_tree_map are dictionaries that map digests to the corresponding WDL object. These are kept as different variables
        because the relationship of objects in/not in the map depend on the context, i.e. caller_tree_map only contains Tasks, Calls, and Workflows that call tasks
        tree_map contains all nodes in the tree
    tree_parser.hierarchy is a hierarchical struct representation of the tree

    Example output stats:
        print("Unique tree paths in WDL:")
        pprint(tree_parser.level.get_all_paths())
        print(f"Counts for parsed node types: {tree_parser.stats}")
        print(f"Number of input nodes in the parsed tree: {len(tree_parser.tree)}")
        print(f"Number of nodes placed into hierarchy: {count_nodes(hierarchy[0])}")
"""

class EncodeSetToList(json.JSONEncoder):
    """Custom encoder to convert sets to lists"""
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return super().default(obj)

def compute_call_digest(call):
    """Calls and Scatters do not have a digest available, compute the MD5 hash of critical keys (id + position information)"""
    call_json = json.dumps([call.workflow_node_id, call.pos], sort_keys=True, cls=EncodeSetToList)
    md5_hash = hashlib.md5(call_json.encode()).hexdigest()
    return md5_hash

def create_hierarchy(nodes):
    # Add the root node (which must be first)
    tree = [nodes[0].copy()]
    assert nodes[0]['parent_digest'] is None, "Root node must be first in the list"

    for node in nodes[1:]:
        node_position = node['parent_digest']
        # node_position = node['metadata_digest']
        # For every entry in node_position list, find the corresponding node
        def place_node_by_position(tree, cur_position, new_node):
            # Is the cur_position on the current level of the tree?
            found_node = [node for node in tree if node['digest'] == cur_position[0]]
            # if len(found_node) > 1:
            #     print("Found node", found_node, cur_position)
            if found_node and len(cur_position) > 1:
                assert len(found_node) == 1, "Found multiple nodes with the same name"
                # Search the children of this node:
                return place_node_by_position(found_node[0]['children'], cur_position[1:], new_node)
            elif found_node and len(cur_position) == 1:
                # We are at the correct node, add the new node as a child
                assert len(found_node) == 1, "Found multiple nodes with the same name"
                if 'children' not in found_node[0]:
                    found_node[0]['children'] = []
                found_node[0]['children'].append(new_node.copy())
                return tree
            else:
                assert False, f"Node parent not found for child node: {new_node}"
        place_node_by_position(tree, node_position, node)
    return tree

def count_nodes(node):
    # Count the current node
    count = 1
    # If the node has children, count them recursively
    if "children" in node:
        for child in node["children"]:
            count += count_nodes(child)
    return count

class WDL_Task_Tree:
    """Recursively walk through the WDL tree and store the results in a tree list"""
    def __init__(self, tree=[]):
        self.tree = tree
        self.hierarchy = []
        self.level = Level()
        self.capture_inputs = True
        self.print_structure = False  # Print task structure + conditionals
        self.tree_map = {}
        self.caller_tree_map = {}
        self.stats = {}  # {'sub-workflow': 5, 'root-workflow': 1, 'task': 39, 'scatter': 5}

    # Graph node objects:
    # WDL.Tree.Calls are wrappers for Tasks or Workflows and are never nodes
    # workflows can be a root workflow (no caller) or sub-workflow (called by a call)
    # Gathers are implicit objects from scatters and are not nodes
    node_objs = [WDL.Tree.Workflow, WDL.Tree.Scatter, WDL.Tree.Call, WDL.Tree.Task]
    stat_types = ['sub-workflow', 'root-workflow', 'task', 'scatter']

    def parse_root_workflow(self, workflow):
        """Main function to start parsing a root workflow (a workflow with no parent)"""
        if not isinstance(workflow, WDL.Tree.Workflow):
            log.error(f"parse_root_workflow() called with non-workflow: '{workflow}'")
            return
        level = self.level
        if self.print_structure:
            log.info(f"{level.get_name()}Root Workflow: {workflow.name}")
        self.add_node(workflow)
        # workflow.body contains all the calls in the workflow
        with level.enter_level(workflow.digest, workflow.name):
            for idx, call in enumerate(workflow.body):
                self.parse_tree_objs(call, level)
        # Add the metadata paths to the tree:
        self.add_metadata_paths()
        # After we finish parsing the workflow into a tree we create a hierarchical struct from 'tree'
        self.hierarchy = create_hierarchy(self.tree)
        self.compute_node_stats()  # Update stats for what we parsed
        return

    def compute_node_stats(self):
        node_types = set([x['type'] for x in self.tree])
        self.stats = {node_type: sum(1 for x in self.tree if x['type'] == node_type) for node_type in node_types}
        # Add 0 counts for any missing keys:
        for key in self.stat_types:
            if key not in self.stats:
                self.stats[key] = 0

    def add_metadata_paths(self):
        """
        Add the metadata path, that is the namespace path we would see in a metadata file or would
        be used to set a variable in an inputs.json file. Such as:
            root-workflow.sub-workflow.(sub-workflow2).variable_name
            WholeGenomeGermlineSingleSample.sample_and_unmapped_bams.sample_name
        (note this does not hide scatter nodes which are not real as far as WDL cares)
        """
        level_digest_paths = self.level.get_all_digest_paths()

        # WDL does not consider scatter objects as namespaces, so we need to remove them from the paths
        metadata_digest_paths = []
        valid_digests = [node['digest'] for node in self.tree if node['type'] != "scatter"]
        for digest_path in level_digest_paths:
            metadata_digest_paths.append([d for d in digest_path if d in valid_digests])

        # Given a set of valid digest paths add a metadata path for each node
        for node in self.tree:
            if not node['parent_digest']:
                node['metadata_digest'] = None
                continue
            digest_index = [i for i, x in enumerate(level_digest_paths) if x == node['parent_digest']]
            assert len(digest_index) == 1, "There should be only a single matching digest"
            node['metadata_digest'] = metadata_digest_paths[digest_index[0]]
        return

    # Recursively retrieve all parent conditional statements until a node_objs is reached
    def retrieve_parent_conditions(self, caller, conditions=None):
        if conditions is None:
            conditions = []
        if type(caller.parent) is WDL.Tree.Conditional:
            conditional = str(caller.parent.expr)
            if type(conditional) is not str:
                conditional = str(conditional)
            assert isinstance(conditional, str), "Conditional is not a string"
            conditions.append(conditional)
            if type(caller.parent.parent) not in WDL_Task_Tree.node_objs: # type: ignore
                conditions = self.retrieve_parent_conditions(caller.parent, conditions)
        return conditions

    def add_node_metadata(self, level, type, alias, name, position, digest, conditions, other_fields=None):
        node = {
            "type": type,
            "alias": alias,
            "name": name,
            "position": position,
            "conditions": conditions,
            "parent": level.get_level() if level else None,
            "digest": digest,
            "parent_digest": level.get_level_as_digest() if level else None
        }
        if other_fields:
            # Add the dictionary in other_fields to the node:
            node.update(other_fields)
        self.tree.append(node)
        return

    # Add a node to the tree, this function is used instead of directly adding node types
    # to make it easy to change the nodes that are captured
    def add_node(self, node, level=None, caller=None):
        if type(node) not in WDL_Task_Tree.node_objs:
            return
        match type(node):
            case WDL.Tree.Scatter:
                self.add_node_scatter(node, level)
            case WDL.Tree.Workflow:
                self.add_node_workflow(node, level, caller)
            case WDL.Tree.Task:
                self.add_node_task(node, level, caller)
            case WDL.Tree.Call:
                if type(node.callee) is WDL.Tree.Task:
                    self.add_node_task(node.callee, level, node)
                elif type(node.callee) is WDL.Tree.Workflow:
                    self.add_node_workflow(node.callee, level, node)
                else:
                    assert True, "WDL.Tree.Call callee not Task or Workflow"
            case _:
                assert True, f"add_node() called with non-node object: {node}"
        return

    def add_node_workflow(self, workflow, level, caller=None):
        if not caller:
            # Root workflows have no caller (and cannot be aliased, 'as'), and are the top level of the tree
            self.add_node_metadata(level, "root-workflow", workflow.name, workflow.name, workflow.pos, workflow.digest, None)
            self.tree_map[workflow.digest] = workflow
        else:
            # Sub-workflows have a caller, and are nested within the tree. Position is the caller since that's got the context
            conditions = self.retrieve_parent_conditions(caller)
            conditions.reverse()  # Match order in the WDL (not closest parent first)
            workflow_digest = workflow.digest + '__' + compute_call_digest(caller)
            self.add_node_metadata(level, "sub-workflow", caller.name, '.'.join(caller.callee_id), caller.pos, workflow_digest, conditions)
            self.tree_map[workflow_digest] = caller
            # Calls eventually need both the callee and the caller:
            self.caller_tree_map[compute_call_digest(caller)] = caller
            self.caller_tree_map[workflow.digest] = workflow
        return

    def get_scatter_variable(self, scatter):
        if type(scatter.expr) in [WDL.Expr.Get, WDL.Expr.Apply]:
            scatter_variable = f"Scatter({scatter.variable} in {str(scatter.expr)})"
            return scatter_variable
        else:
            assert True, "Scatter expression is not a WDL.Expr.Get"

    def add_node_scatter(self, scatter, level):
        conditions = self.retrieve_parent_conditions(scatter)
        conditions.reverse()  # Match order in the WDL (not closest parent first)
        scatter_variable = self.get_scatter_variable(scatter)
        # Scatter_ + scatter.variable causes name collisions, so we overspecify both name/alias
        scatter_digest = compute_call_digest(scatter)
        self.add_node_metadata(level, "scatter", scatter_variable, scatter_variable, scatter.pos, scatter_digest, conditions)
        self.tree_map[scatter_digest] = scatter
        return

    def add_node_task(self, task, level, caller):
        conditions = self.retrieve_parent_conditions(caller)
        conditions.reverse()  # Match order in the WDL (not closest parent first)
        other_fields = self.return_node_call_inputs(task, caller)  # Add the inputs/inputs_caller fields
        if not other_fields:
            other_fields = {}
        other_fields['callee_name'] = task.name
        other_fields['callee_pos'] = task.pos
        other_fields['callee_digest'] = task.digest
        task_digest = task.digest + '__' + compute_call_digest(caller)
        self.add_node_metadata(level, "task", caller.name, '.'.join(caller.callee_id), caller.pos, task_digest, conditions, other_fields)
        self.tree_map[task_digest] = caller
        # Calls eventually need both the callee and the caller:
        self.caller_tree_map[compute_call_digest(caller)] = caller
        self.caller_tree_map[task.digest] = task
        return

    def return_node_call_inputs(self, task, caller):
        if not self.capture_inputs:
            return None
        inputs = {"inputs_caller": [input for input in caller.inputs]}
        # inputs['inputs_task'] = [input.name for input in task.inputs]
        inputs['inputs'] = [input.name for input in task.available_inputs]
        return inputs

    def return_node_inputs(self, node):
        if self.capture_inputs:
            return {"inputs": [input for input in node.inputs]}
        else:
            return None

    def parse_tree_objs(self, call, level):
        """Parse functions handle printing and descending through the tree, add functions do the work
        of adding the nodes to the tree structure"""
        match type(call):
            case WDL.Tree.Call:
                if type(call.callee) is WDL.Tree.Task:
                    self.parse_task(call.callee, level, call)
                elif type(call.callee) is WDL.Tree.Workflow:
                    self.parse_workflow(call.callee, level, call)
                else:
                    assert True, "WDL.Tree.Call callee not Task or Workflow"
            case WDL.Tree.Conditional:
                self.parse_conditional(call, level)
            case WDL.Tree.Scatter:
                if self.print_structure:
                    log.info(f"{level.get_name()} Scatter: {call.variable} {str(call.expr)}")
                self.parse_scatter(call, level)
            case WDL.Tree.Workflow:
                if self.print_structure:
                    log.info(f"{level.get_name()} Workflow: {call.name}")
                self.parse_workflow(call, level)
            case WDL.Tree.Decl:
                pass  # Ignore strings and other objects
            case WDL.Tree.Task:
                if self.print_structure:
                    log.info(f"{level.get_name()} Task: {call.name}")
                assert True, "Unhanded task"
            case _:
                if self.print_structure:
                    log.info(f"{level.get_name()} Other: {call.type} {call.name}")

    def parse_workflow(self, workflow, level, caller=None):
        self.add_node(workflow, level, caller)
        # If we have a caller we should use the caller name otherwise use the workflow name
        workflow_name = workflow.name if not caller else caller.name
        workflow_digest = workflow.digest if not caller else workflow.digest + '__' + compute_call_digest(caller)
        with level.enter_level(workflow_digest, workflow_name):
            for call in workflow.body:
                self.parse_tree_objs(call, level)
        return

    def parse_conditional(self, call, level):
        for body in call.body:
            self.parse_tree_objs(body, level)
        return

    def parse_scatter(self, scatter, level):
        self.add_node(scatter, level)

        scatter_variable = self.get_scatter_variable(scatter)
        with level.enter_level(compute_call_digest(scatter), scatter_variable):
            for body in scatter.body:
                self.parse_tree_objs(body, level)
        return

    def parse_task(self, task, level, caller):
        self.add_node(task, level, caller)
        conditions = self.retrieve_parent_conditions(caller)
        conditions.reverse()  # Match order in the WDL (not closest parent first)
        if conditions and self.print_structure:
            log.info(f"{level.get_name()}       Conditions: {conditions}")
        if self.print_structure:
            log.info(f"{level.get_name()} Task: {caller.name} -> {task.name} [{'.'.join(caller.callee_id)}]")
        return

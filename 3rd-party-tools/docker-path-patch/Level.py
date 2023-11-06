from contextlib import contextmanager

class Level:
    """
    contextmanager to track the current recursive level. We use the digest (or workflow_node_id) to track
    the current recursion depth into the WDL AST tree (to avoid collisions). And provide a way to get back
    to a human readable name for the current level.

    Args:
        name: human readable name of the current level
        digest: WDL digest for the current node or workflow_node_id for a sub-workflow

    Attributes:
        all_paths is a list of the mapped levels encountered during recursion
    """
    def __init__(self, name=None):
        self.name = name
        self.digests = []
        self.mapping = {}
        self.current_level = []
        self.all_paths = []
        self.all_digest_paths = []

    def add_digest(self, digest, name):
        self.digests.append(digest)
        self.mapping[digest] = name
        self.current_level.append(digest)

    @contextmanager
    def enter_level(self, digest, name):
        self.add_digest(digest, name)
        self.all_paths.append(self.get_name())
        self.all_digest_paths.append(self.get_level_as_digest())
        try:
            yield
        finally:
            self.current_level.pop()

    def get_level(self):
        # Returns the current level as a list of mapped digests
        return [self.mapping[digest] for digest in self.current_level]

    def get_level_as_digest(self):
        # Returns the current level as a list of mapped digests
        return self.current_level.copy()

    def get_name(self):
        # Returns a string representation of the current level separated by '.'
        return '.'.join([self.mapping[digest] for digest in self.current_level])

    def get_all_paths(self):
        return self.all_paths

    def get_all_digest_paths(self):
        return self.all_digest_paths

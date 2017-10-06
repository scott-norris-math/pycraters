class Parameters(dict):

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)


    def __getitem__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)


    def __setattr__(self, name, value):
        if name in self:
            self[name] = value
        else:
            raise AttributeError("No such attribute: " + name)


    def __setitem__(self, name, value):
        if name in self:
            self[name] = value
        else:
            raise AttributeError("No such attribute: " + name)


    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)

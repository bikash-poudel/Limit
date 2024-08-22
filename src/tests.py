
class flux_connection:
    def __init__(self, left_node, right_node):
        self.left_node = left_node
        self.right_node = right_node

        # register connection
        self.left_node.RegisterConnection(newConnection=self)
        self.right_node.RegisterConnection(newConnection=self)


class LiquidDiffusion(flux_connection):


    def method1(self):
        return "Method 1 of LiquidDiffusion"

    def method2(self):
        return "Method 2 of LiquidDiffusion"




class VaporDiffusion(flux_connection):
    def method1(self):
        return "Method 1 of VaporDiffusion"

    def method2(self):
        return "Method 2 of VaporDiffusion"


class flux_node(object):
    Dict_of_all_flux_nodes = {}
    current_ID = 1

    def __init__(self, conc_iso_liquid={"2H": 1.0, "18O": 1.0}, T=283.15):
        self.__ID = self.__class__.current_ID
        self.__class__.current_ID += 1
        self.__class__.Dict_of_all_flux_nodes[self.__ID] = self

        self.__Connections = []
        self.__Connections_to_iso_storages = []
        self.__Connections_to_boundaries = []
        self.__conc_iso_liquid = conc_iso_liquid
        self.T = T

    def __getattr__(self, name):
        for connection in self.__Connections:
            if connection.__class__.__name__ == name:
                setattr(self, name, connection)  # Dynamically attach the connection instance
                return connection
        raise AttributeError(f"No connection named {name}")

    def get_connections(self):
        return self.__Connections

    connections = property(get_connections, None, None, "List of all flux connections connected with this node")

    def RegisterConnection(self, newConnection):
        """
        Registers the given connection.
        """
        self.__Connections.append(newConnection)

        if newConnection.left_node == self:
            other_node = newConnection.right_node
        elif newConnection.right_node == self:
            other_node = newConnection.left_node
        else:
            raise NotImplementedError

        if isinstance(other_node, flux_node):
            other_node.__Connections.append(newConnection)
        elif isinstance(other_node, iso_storage):
            self.__Connections_to_iso_storages.append(other_node)
        elif isinstance(other_node, iso_atmosphere):
            self.__Connections_to_boundaries.append(other_node)


# Example classes to use with flux_node for testing
class iso_storage:
    pass

class iso_atmosphere:
    pass

# Create flux nodes
node1 = flux_node()
node2 = flux_node()

# Create instances of connection classes
ld_instance = LiquidDiffusion(node1, node2)
vd_instance = VaporDiffusion(node1, node2)

# Accessing methods through dynamic attribute access
print(node1.LiquidDiffusion.method1())  # Output: Method 1 of LiquidDiffusion
print(node1.VaporDiffusion.method2())   # Output: Method 2 of VaporDiffusion


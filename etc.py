class SoilLayer:
    """Represents a single soil layer with position and depth"""

    def __init__(self, layer_id, top_depth, bottom_depth, properties=None):
        self.layer_id = layer_id
        self.top_depth = top_depth
        self.bottom_depth = bottom_depth
        self.thickness = bottom_depth - top_depth
        self.properties = properties or {}
        self.connections = set()  # Store connections to other layers

    def __repr__(self):
        return f"Layer({self.layer_id}, {self.top_depth}-{self.bottom_depth}m)"

    def connects_to(self, other_layer):
        """Check if this layer connects to another layer based on depth overlap"""
        return not (self.bottom_depth <= other_layer.top_depth or
                    self.top_depth >= other_layer.bottom_depth)


class SoilCell:
    """Represents a horizontal cell containing vertical layers"""

    def __init__(self, cell_id, x, y):
        self.cell_id = cell_id
        self.x = x
        self.y = y
        self.layers = []  # List of SoilLayer objects, ordered by depth

    def add_layer(self, layer):
        """Add a layer to the cell, maintaining depth order"""
        self.layers.append(layer)
        # Sort layers by depth
        self.layers.sort(key=lambda l: l.top_depth)

    def get_layer_at_depth(self, depth):
        """Get the layer that contains the given depth"""
        for layer in self.layers:
            if layer.top_depth <= depth < layer.bottom_depth:
                return layer
        return None

    def __repr__(self):
        return f"Cell({self.cell_id}, [{self.x},{self.y}], {len(self.layers)} layers)"


class SoilDomain:
    """Manages the entire soil domain with cells and layer connections"""

    def __init__(self):
        self.cells = {}  # Dict with (x,y) as key
        self.connections = []  # List of all connections

    def add_cell(self, cell):
        """Add a cell to the domain"""
        self.cells[(cell.x, cell.y)] = cell

    def get_adjacent_cells(self, cell):
        """Get adjacent cells (up, down, left, right)"""
        adjacent = []
        directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # N, E, S, W

        for dx, dy in directions:
            neighbor_pos = (cell.x + dx, cell.y + dy)
            if neighbor_pos in self.cells:
                adjacent.append(self.cells[neighbor_pos])

        return adjacent

    def connect_layers_between_cells(self, cell1, cell2):
        """Connect layers between two adjacent cells"""
        connections = []

        # For each layer in cell1, find overlapping layers in cell2
        for layer1 in cell1.layers:
            for layer2 in cell2.layers:
                if layer1.connects_to(layer2):
                    # Calculate overlap percentage
                    overlap_top = max(layer1.top_depth, layer2.top_depth)
                    overlap_bottom = min(layer1.bottom_depth, layer2.bottom_depth)
                    overlap_thickness = overlap_bottom - overlap_top

                    if overlap_thickness > 0:
                        connection = {
                            'cell1': cell1.cell_id,
                            'cell2': cell2.cell_id,
                            'layer1': layer1.layer_id,
                            'layer2': layer2.layer_id,
                            'overlap_top': overlap_top,
                            'overlap_bottom': overlap_bottom,
                            'overlap_thickness': overlap_thickness,
                            'connection_strength': min(
                                overlap_thickness / layer1.thickness,
                                overlap_thickness / layer2.thickness
                            )
                        }
                        connections.append(connection)

                        # Store connection in layer objects
                        layer1.connections.add((cell2.cell_id, layer2.layer_id))
                        layer2.connections.add((cell1.cell_id, layer1.layer_id))

        return connections

    def build_all_connections(self):
        """Build connections between all adjacent cells in the domain"""
        self.connections = []

        for cell in self.cells.values():
            adjacent_cells = self.get_adjacent_cells(cell)
            for adjacent_cell in adjacent_cells:
                # Only connect each pair once
                if cell.cell_id < adjacent_cell.cell_id:
                    connections = self.connect_layers_between_cells(cell, adjacent_cell)
                    self.connections.extend(connections)

        return self.connections

    def visualize_connections(self, cell_id=None):
        """Visualize connections for a specific cell or all cells"""
        if cell_id:
            # Show connections for specific cell
            print(f"\nConnections for Cell {cell_id}:")
            for conn in self.connections:
                if conn['cell1'] == cell_id or conn['cell2'] == cell_id:
                    other_cell = conn['cell2'] if conn['cell1'] == cell_id else conn['cell1']
                    other_layer = conn['layer2'] if conn['cell1'] == cell_id else conn['layer1']
                    print(f"  → Cell {other_cell}, Layer {other_layer}: "
                          f"Overlap {conn['overlap_top']}-{conn['overlap_bottom']}m "
                          f"(strength: {conn['connection_strength']:.2f})")
        else:
            # Show all connections
            print("\nAll connections in domain:")
            for conn in self.connections:
                print(f"Cell {conn['cell1']} Layer {conn['layer1']} ↔ "
                      f"Cell {conn['cell2']} Layer {conn['layer2']}: "
                      f"{conn['overlap_top']}-{conn['overlap_bottom']}m")


# Example usage and demonstration
def create_example_domain():
    """Create an example soil domain for testing"""
    domain = SoilDomain()

    # Create cells
    cell1 = SoilCell("C1", 0, 0)
    cell2 = SoilCell("C2", 1, 0)
    cell3 = SoilCell("C3", 0, 1)

    # Add layers to cell1 (uniform layers)
    cell1.add_layer(SoilLayer("L1", 0, 2, {"type": "topsoil"}))
    cell1.add_layer(SoilLayer("L2", 2, 5, {"type": "clay"}))
    cell1.add_layer(SoilLayer("L3", 5, 8, {"type": "sand"}))

    # Add layers to cell2 (different depths)
    cell2.add_layer(SoilLayer("L1", 0, 1.5, {"type": "topsoil"}))
    cell2.add_layer(SoilLayer("L2", 1.5, 4, {"type": "clay"}))
    cell2.add_layer(SoilLayer("L3", 4, 7, {"type": "sand"}))
    cell2.add_layer(SoilLayer("L4", 7, 9, {"type": "gravel"}))

    # Add layers to cell3 (yet another configuration)
    cell3.add_layer(SoilLayer("L1", 0, 3, {"type": "topsoil"}))
    cell3.add_layer(SoilLayer("L2", 3, 6, {"type": "clay"}))
    cell3.add_layer(SoilLayer("L3", 6, 10, {"type": "bedrock"}))

    # Add cells to domain
    domain.add_cell(cell1)
    domain.add_cell(cell2)
    domain.add_cell(cell3)

    return domain


# Demonstration
if __name__ == "__main__":
    # Create example domain
    domain = create_example_domain()

    # Build all connections
    connections = domain.build_all_connections()

    print("=== Soil Domain Layer Connection System ===")
    print(f"Total cells: {len(domain.cells)}")
    print(f"Total connections found: {len(connections)}")

    # Show domain structure
    print("\n=== Domain Structure ===")
    for cell in domain.cells.values():
        print(f"\n{cell}:")
        for layer in cell.layers:
            print(f"  {layer}")

    # Show connections
    domain.visualize_connections()

    # Show specific cell connections
    domain.visualize_connections("C1")
    domain.visualize_connections("C2")

    # Example: Find which layers connect at a specific depth
    print("\n=== Connection Analysis at Depth 2m ===")
    depth = 2.0
    for cell in domain.cells.values():
        layer = cell.get_layer_at_depth(depth)
        if layer:
            print(f"Cell {cell.cell_id}: {layer} contains depth {depth}m")
            if layer.connections:
                print(f"  Connected to: {layer.connections}")
import pysmile
import pysmile_license


# Helper to create a CPT node in the network
def create_cpt_node(net, id, name, outcomes, x_pos, y_pos):
    handle = net.add_node(pysmile.NodeType.CPT, id)
    net.set_node_name(handle, name)
    net.set_node_position(handle, x_pos, y_pos, 85, 55)

    initial_outcome_count = net.get_outcome_count(handle)

    # Fill/rename existing default outcomes
    for i in range(0, min(initial_outcome_count, len(outcomes))):
        net.set_outcome_id(handle, i, outcomes[i])

    # If more outcomes are required, add them
    for i in range(initial_outcome_count, len(outcomes)):
        net.add_outcome(handle, outcomes[i])

    return handle

def build_structure_only(save_path="BN.xdsl"):
    net = pysmile.Network()

    # ---- Nodes & states (no probabilities) ----
    # PGA: 20 bins PGA_00..PGA_19
    pga_states = [f"PGA_{i:02d}" for i in range(20)]
    create_cpt_node(net, "PGA", "Peak Ground Acceleration", pga_states, 160, 40)

    # AGE: 6 levels
    age_states = [f"AGE_{a}" for a in (0, 10, 20, 30, 40, 50)]
    create_cpt_node(net, "AGE", "Reactor Age", age_states, 60, 40)

    # LOCA: binary
    create_cpt_node(net, "LOCA", "Loss of Coolant Accident", ["NO_LOCA", "YES_LOCA"], 110, 140)

    # ACS1: binary
    create_cpt_node(net, "ACS1", "Auxiliary Cooling System 1 Failure", ["NO_ACS1", "YES_ACS1"], 360, 140)

    # Discretized time bins (labels only; pick what you like)
    edges = [0, 1, 5, 30, 120, 300, 600, 900, 1200, float("inf")]
    tbin_labels = [f"[{int(edges[k])},{'inf' if edges[k+1]==float('inf') else int(edges[k+1])})"
                   for k in range(len(edges)-1)]

    create_cpt_node(net, "t_loca_bin", "Time to LOCA (binned)", tbin_labels, 110, 240)
    create_cpt_node(net, "t_acs1_bin", "Time to ACS1 fail (binned)", tbin_labels, 360, 240)

    # ---- Structure only (arcs), no CPTs ----
    net.add_arc("PGA", "LOCA")
    net.add_arc("AGE", "LOCA")

    net.add_arc("PGA", "ACS1")
    net.add_arc("AGE", "ACS1")

    net.add_arc("LOCA", "t_loca_bin")
    net.add_arc("ACS1", "t_acs1_bin")

    # Save the empty-parameter structure (you can fill CPTs later in GeNIe or via code)
    net.write_file(save_path)
    print(f"Structure-only network saved to {save_path}")

if __name__ == "__main__":
    build_structure_only()
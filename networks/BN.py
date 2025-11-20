import pysmile
import pysmile_license


def create_cpt_node(net, node_id, name, outcomes, x_pos, y_pos, width=85, height=55, registry=None):
    """
    Utility helper that creates a CPT node (structure only) with the desired set of outcomes.
    """
    handle = net.add_node(pysmile.NodeType.CPT, node_id)
    net.set_node_name(handle, name)
    net.set_node_position(handle, x_pos, y_pos, width, height)

    existing = net.get_outcome_count(handle)
    for idx in range(min(existing, len(outcomes))):
        net.set_outcome_id(handle, idx, outcomes[idx])

    for idx in range(existing, len(outcomes)):
        net.add_outcome(handle, outcomes[idx])

    if registry is not None:
        registry.add(node_id)

    return handle


def time_bin_labels(edges):
    labels = []
    for start, stop in zip(edges[:-1], edges[1:]):
        right = "inf" if stop == float("inf") else int(stop)
        labels.append(f"[{int(start)},{right})")
    return labels


def build_structure_only(save_path="BN.xdsl"):
    net = pysmile.Network()
    node_registry = set()

    # --- state definitions -------------------------------------------------
    pga_states = [f"PGA_{i:02d}" for i in range(20)]
    age_states = [f"AGE_{i}0" for i in range(6)]  # AGE_00 .. AGE_50
    # Coarser bins keep CPT sizes manageable while preserving the causal structure.
    event_time_bins = time_bin_labels([0, 30, 300, float("inf")])
    repair_time_bins = time_bin_labels([0, 60, 240, float("inf")])
    distance_states = ["Distance_250", "Distance_500", "Distance_750", "Distance_1000"]

    def yes_no_states(base):
        return [f"NO_{base}", f"YES_{base}"]

    # --- node placement map (rough layout to keep the diagram readable) ----
    positions = {
        "AGE": (80, 40), "PGA": (220, 40), "DISTANCE": (520, 40),
        "LOCA": (150, 140), "t_loca": (150, 240),
        "LOOP": (260, 140), "t_loop": (260, 240),
        "MSLB": (370, 140), "t_mslb": (370, 240),
        "MSLBH2": (480, 140), "t_MSLBH2": (480, 240),
        "EXP": (590, 140),
        "OC": (330, 40),
        "LOOPH2": (700, 140), "t_LOOPH2": (700, 240),
        "LHS": (810, 140), "t_LHS": (810, 240),
        "ACS": (40, 140), "t_acs": (40, 240), "rt_acs": (40, 340),
        "EDG": (150, 340), "t_edg": (150, 440), "rt_edg": (150, 540),
        "PDP": (260, 340), "t_pdp": (260, 440), "rt_pdp": (260, 540),
        "Reactor": (520, 440),
    }

    # --- discrete equipment / scenario nodes -------------------------------
    create_cpt_node(net, "PGA", "PGA", pga_states, *positions["PGA"], registry=node_registry)
    create_cpt_node(net, "AGE", "AGE", age_states, *positions["AGE"], registry=node_registry)
    create_cpt_node(net, "DISTANCE", "DISTANCE", distance_states, *positions["DISTANCE"], registry=node_registry)

    binary_nodes = [
        "LOCA", "LOOP", "MSLB", "MSLBH2", "EXP", "OC", "LOOPH2", "LHS", "ACS", "EDG", "PDP",
    ]

    for node_id in binary_nodes:
        create_cpt_node(
            net,
            node_id,
            node_id,
            yes_no_states(node_id),
            *positions.get(node_id, (0, 0)),
            registry=node_registry,
        )

    # --- time nodes --------------------------------------------------------
    time_nodes = [
        "t_loca",
        "t_loop",
        "t_mslb",
        "t_MSLBH2",
        "t_LOOPH2",
        "t_LHS",
        "t_acs",
        "t_edg",
        "t_pdp",
    ]

    for node_id in time_nodes:
        create_cpt_node(
            net,
            node_id,
            node_id,
            event_time_bins,
            *positions.get(node_id, (0, 0)),
            registry=node_registry,
        )

    residual_time_nodes = [
        "rt_acs",
        "rt_edg",
        "rt_pdp",
    ]

    for node_id in residual_time_nodes:
        create_cpt_node(
            net,
            node_id,
            node_id,
            repair_time_bins,
            *positions.get(node_id, (0, 0)),
            registry=node_registry,
        )

    # Reactor performance node (placeholder states for structure only)
    create_cpt_node(
        net,
        "Reactor",
        "Reactor",
        ["SAFE", "UNSAFE"],
        *positions["Reactor"],
        registry=node_registry,
    )

    # --- network structure -------------------------------------------------
    arcs = [
        ("AGE", "LOCA"), ("PGA", "LOCA"), ("LOCA", "t_loca"),
        ("PGA", "LOOP"), ("LOOP", "t_loop"),
        ("AGE", "MSLB"), ("PGA", "MSLB"), ("MSLB", "t_mslb"),
        ("PGA", "MSLBH2"), ("MSLBH2", "t_MSLBH2"),
        ("MSLB", "EXP"), ("PGA", "EXP"),
        ("EXP", "LOOPH2"), ("LOOPH2", "t_LOOPH2"),
        ("PGA", "OC"),
        ("EXP", "LHS"), ("OC", "LHS"), ("LHS", "t_LHS"),
        ("AGE", "ACS"), ("PGA", "ACS"), ("ACS", "t_acs"),
        ("AGE", "EDG"), ("PGA", "EDG"), ("EDG", "t_edg"),
        ("AGE", "PDP"), ("PGA", "PDP"), ("PDP", "t_pdp"),
        ("DISTANCE", "LOOPH2"), ("DISTANCE", "LHS"),
        ("t_loca", "Reactor"), ("t_loop", "Reactor"), ("t_mslb", "Reactor"),
        ("t_MSLBH2", "Reactor"), ("t_LOOPH2", "Reactor"), ("t_LHS", "Reactor"),
        ("t_acs", "Reactor"), ("rt_acs", "Reactor"),
        ("t_edg", "Reactor"), ("rt_edg", "Reactor"),
        ("t_pdp", "Reactor"), ("rt_pdp", "Reactor"),
    ]

    for parent, child in arcs:
        if parent not in node_registry or child not in node_registry:
            raise ValueError(f"Arc {parent}->{child} references undefined node(s)")
        try:
            net.add_arc(parent, child)
        except pysmile.SMILEException as exc:
            raise RuntimeError(f"Failed to add arc {parent}->{child}") from exc

    net.write_file(save_path)
    print(f"Structure-only network saved to {save_path}")


if __name__ == "__main__":
    build_structure_only()

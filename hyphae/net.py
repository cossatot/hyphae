from collections import defaultdict
from itertools import product, combinations, chain


import numpy as np
import pandas as pd


def get_connected_subgraphs(graph):
    def backtrack(start, current_set):
        result.append(current_set.copy())

        for vertex in graph:
            if vertex not in current_set:
                if any(neighbor in current_set for neighbor in graph[vertex]):
                    current_set.add(vertex)
                    backtrack(vertex, current_set)
                    current_set.remove(vertex)

    result = []
    for start_vertex in graph:
        backtrack(start_vertex, {start_vertex})

    # convert to list, remove duplicates and sort for consistent output
    unique_subgraphs = list(map(frozenset, result))
    unique_subgraphs = sorted(
        set(unique_subgraphs), key=lambda x: (len(x), tuple(sorted(x)))
    )
    unique_subgraphs = [list(sg) for sg in unique_subgraphs]
    return unique_subgraphs


def pdf_overlap(*eqs):
    overlap = pd.concat(eqs, axis=1).dropna().min(axis=1)
    overlap = overlap[overlap > 0.0]
    return overlap


def intersect_intervals(intervals1, intervals2):
    result = []
    for start1, end1 in intervals1:
        for start2, end2 in intervals2:
            intersection_start = min(start1, start2)
            intersection_end = max(end1, end2)
            if intersection_start >= intersection_end:
                result.append((intersection_start, intersection_end))

    return result


def get_eq_combinations(trench_list, eq_trench_dates):
    eq_lists = [list(eq_trench_dates[trench].keys()) for trench in trench_list]
    return product(*eq_lists)


def find_common_intervals(subgraphs, eq_trench_dates):
    common_intervals = {}

    for subgraph in subgraphs:
        if len(subgraph) == 1:
            eq_ages = {
                (k,): v for k, v in eq_trench_dates[subgraph[0]].items()
            }
            common_intervals[tuple(subgraph)] = eq_ages
        else:
            subgraph_results = {}
            eq_combinations = get_eq_combinations(subgraph, eq_trench_dates)

            for combination in eq_combinations:
                pdfs = [
                    eq_trench_dates[trench][eq]
                    for trench, eq in zip(subgraph, combination)
                ]
                overlap = pdf_overlap(*pdfs)
                if len(overlap) > 0:
                    combo_key = tuple(
                        sorted(
                            [eq for trench, eq in zip(subgraph, combination)]
                        )
                    )
                    subgraph_results[combo_key] = overlap
            if len(subgraph_results) > 0:
                common_intervals[tuple(sorted(subgraph))] = subgraph_results

    return common_intervals


def pdf_mean(pdf):
    pdf_mean = np.sum(
        [np.product(age_prob) for age_prob in zip(pdf.index, pdf.values)]
    )
    pdf_mean = pdf_mean / np.trapz(pdf.values, x=pdf.index)
    return pdf_mean


def overlap_pdf_integral(pdf):
    return np.trapz(pdf.values, pdf.index)


def haversine_distance(
    lon_0: float = None,
    lat_0: float = None,
    lon_1: float = None,
    lat_1: float = None,
    R: float = 6371,
) -> float:
    """
    Calculates the great circle distance between two points in lon, lat
    using the haversine formula.
    """
    r_lon_0, r_lon_1, r_lat_0, r_lat_1 = np.radians(
        (lon_0, lon_1, lat_0, lat_1)
    )
    term_1 = np.sin((r_lat_1 - r_lat_0) / 2.0) ** 2
    term_2 = np.cos(r_lat_0) * np.cos(r_lat_1)
    term_3 = np.sin((r_lon_1 - r_lon_0) / 2.0) ** 2

    return 2 * R * np.arcsin(np.sqrt(term_1 + term_2 * term_3))


def get_adjacent_distances(graph, sites):
    """
    Given a distance-weighted graph and a set of sites, returns the distances
    between all pairs of adjacent sites in the set.

    Args:
        graph (dict): Distance-weighted graph where graph[site1][site2] gives
        distance sites (tuple/list/set): Collection of sites to check for
        adjacency

    Returns:
        adjacent_distances: dict of (site1, site2): distance items for all
        adjacent pairs
    """
    adjacent_distances = {}
    sites = set(sites)  # Convert to set for efficient lookup

    for site1 in sites:
        for site2 in graph[site1]:
            if (site2, site1) not in adjacent_distances:
                adjacent_distances[(site1, site2)] = graph[site1][site2]

    return adjacent_distances


def get_displacement_diffs(
    event_sites, event_displacements, displacement_lookup, event_dists
):
    displacement_diffs = {}

    event_connection_graph = {}
    for site1, site2 in event_dists.keys():
        if site1 not in event_connection_graph:
            event_connection_graph[site1] = [site2]
        else:
            if site2 not in event_connection_graph[site1]:
                event_connection_graph[site1].append(site2)

        if site2 not in event_connection_graph:
            event_connection_graph[site2] = [site1]
        else:
            if site1 not in event_connection_graph[site2]:
                event_connection_graph[site2].append(site1)

    for event_site_1 in event_sites:
        site_1 = displacement_lookup[event_site_1]
        connected_sites = set(event_connection_graph[site_1])
        done_sites = set()
        for event_site_2 in event_sites:
            site_2 = displacement_lookup[event_site_2]

            if site_1 == site_2:
                pass
            elif (event_site_2, event_site_1) in displacement_diffs:
                # this pair already done
                connected_sites.remove(site_2)
            elif site_2 not in connected_sites:
                pass  # not connected
            elif site_2 in connected_sites:
                displacement_diffs[(event_site_1, event_site_2)] = np.abs(
                    event_displacements[site_1][event_site_1]
                    - event_displacements[site_2][event_site_2]
                )
                connected_sites.remove(site_2)
        # now do the connected sites w/o displacement
        for site in connected_sites:
            displacement_diffs[(event_site_1, site)] = event_displacements[
                site_1
            ][event_site_1]

    return displacement_diffs


def get_event_displacement_gradients(
    event_sites, displacements, distance_graph
):

    displacement_lookup = {
        event: site
        for site, events in displacements.items()
        for event in events
    }

    sites = [displacement_lookup[site] for site in event_sites]
    event_dists = get_adjacent_distances(distance_graph, sites)

    displacement_diffs = get_displacement_diffs(
        event_sites, displacements, displacement_lookup, event_dists
    )

    displacement_gradients = {}
    for pair, disp_diff in displacement_diffs.items():
        disp_diff_km = 0.001 * disp_diff

        site_1 = displacement_lookup.get(pair[0], pair[0])
        site_2 = displacement_lookup.get(pair[1], pair[1])

        site_pair = (site_1, site_2)
        try:
            site_dist = event_dists[site_pair]
        except KeyError:
            site_dist = event_dists[site_pair[::-1]]

        displacement_gradients[pair] = disp_diff_km / site_dist

    return displacement_gradients


def generate_ordered_sequences(trench_events, overlaps):
    def is_valid_group(group):
        if len(group) == 1:
            return True
        return tuple(sorted(group)) in overlaps or any(
            all(event in overlap for event in group) for overlap in overlaps
        )

    def recursive_group(current_indices, current_sequence):
        # if all events are used, it's a valid sequence
        if all(
            i == len(events)
            for i, events in zip(current_indices, trench_events.values())
        ):
            yield current_sequence
            return

        # try grouping the event from each trench
        for trench_index, (trench_events) in enumerate(trench_events.items()):
            if current_indices[trench_index] < len(trench_events):
                new_group = [trench_events[current_indices[trench_index]]]
                new_indices = list(current_indices)
                new_indices[trench_index] += 1

                # try adding events from other trenches to this group
                for other_index, (other_trench, other_events) in enumerate(
                    trench_events.items()
                ):
                    if other_index != trench_index and (
                        current_indices[other_index] < len(other_events)
                    ):
                        potential_group = new_group + [
                            other_events[current_indices[other_index]]
                        ]
                        if is_valid_group(potential_group):
                            new_group = potential_group
                            new_indices[other_index] += 1
                        else:
                            break

                # only proceed if the group is valid
                if is_valid_group(new_group):
                    # recursive call iwth the new group added to the sequence
                    yield from recursive_group(
                        tuple(new_indices),
                        current_sequence + [tuple(new_group)],
                    )

    # start with all indices at 0
    initial_indices = tuple(0 for _ in trench_events)
    yield from recursive_group(initial_indices, [])


def generate_earthquake_sequences(trench_events, overlaps):
    def build_graph():
        graph = defaultdict(set)
        for overlap_events in overlaps.keys():
            for event1, event2 in combinations(overlap_events, 2):
                graph[event1].add(event2)
                graph[event2].add(event1)
        return graph

    def is_valid_group(group):
        if len(group) == 1:
            return True
        sorted_group = tuple(sorted(group))
        return sorted_group in overlaps or any(
            set(sorted_group).issubset(set(k)) for k in overlaps.keys()
        )

    def get_trench(event):
        return next(
            trench
            for trench, events in trench_events.items()
            if event in events
        )

    def get_event_index(event):
        trench = get_trench(event)
        return trench_events[trench].index(event)

    graph = build_graph()
    events = [
        event
        for trench_events in trench_events.values()
        for event in trench_events
    ]
    n = len(events)

    def dfs_sequences(index=0, current_sequence=None, used_events=None):
        if current_sequence is None:
            current_sequence = []
        if used_events is None:
            used_events = set()

        if len(used_events) == n:
            yield current_sequence
            return

        current_event = events[index]
        if current_event in used_events:
            yield from dfs_sequences(index + 1, current_sequence, used_events)
            return

        # Try adding to an existing group
        for i, group in enumerate(current_sequence):
            if any(e in graph[current_event] for e in group):
                new_group = group + [current_event]
                if is_valid_group(new_group):
                    # Check stratigraphic order
                    if all(
                        get_event_index(e) < get_event_index(current_event)
                        for e in group
                        if get_trench(e) == get_trench(current_event)
                    ):
                        new_sequence = (
                            current_sequence[:i]
                            + [new_group]
                            + current_sequence[i + 1 :]
                        )
                        new_used = used_events | {current_event}
                        yield from dfs_sequences(
                            index + 1, new_sequence, new_used
                        )

        # Try creating a new group
        new_sequence = current_sequence + [[current_event]]
        new_used = used_events | {current_event}
        yield from dfs_sequences(index + 1, new_sequence, new_used)

    return list(dfs_sequences())

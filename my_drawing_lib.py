from textwrap import wrap
import numpy as np
import pandas as pd
import networkx as nx
import plotly.express as px
import matplotlib.pyplot as plt


def draw_synonyms(drug_bank_id: str, synonyms: pd.DataFrame) -> None:
    syns = synonyms.loc[drug_bank_id, 'synonyms']
    edge_list = [(drug_bank_id, syn) for syn in syns]
    g = nx.from_edgelist(edge_list)

    pos = nx.spring_layout(g, k=1, seed=42)

    # node sizes
    central_size = len(drug_bank_id) * 800
    node_sizes = [central_size if node == drug_bank_id else central_size/2 for node in g]

    # drawing nodes and edges
    plt.figure(figsize=(12, 10))  # Larger area for better readability
    nx.draw_networkx_edges(g, pos, edge_color="gray", width=1.2)
    nx.draw_networkx_nodes(g, pos, node_size=node_sizes, node_color="lightblue", edgecolors="darkblue")

    # labeling the central node (in the middle of the node):
    nx.draw_networkx_labels(
        g,
        {drug_bank_id: pos[drug_bank_id]},  # Only the central node
        {drug_bank_id: drug_bank_id},
        font_size=10,
        verticalalignment="center"
    )

    # labeling other nodes (under the nodes):

    # first we shift the position of other nodes down by 0.2 of an inch
    shifted_pos = {node: (x, y-0.2) for node, (x, y) in pos.items() if node != drug_bank_id}
    nx.draw_networkx_labels(
        g,
        shifted_pos,
        labels={node: node for node in shifted_pos},
        font_size=10,
        verticalalignment="top"
    )

    # showing the graph
    plt.axis("off")
    plt.margins(0.2)
    plt.show()


def draw_bipartite_graph(pathways: pd.DataFrame) -> None:
    # Empty graph
    b = nx.Graph()

    # first we designate the id nodes as left
    left_nodes = list(pathways['smpdb-id'])

    # we add left nodes with the 'bipartite' tag set to 0 to represent one side of the graph
    b.add_nodes_from(left_nodes, bipartite=0)

    # we add the edges, which adds missing right nodes
    edges = [(name, drug) for name, drugs in zip(pathways['smpdb-id'], pathways['drugs']) for drug in drugs]
    b.add_edges_from(edges)

    # calculate a bipartite layout for the nodes
    pos = nx.bipartite_layout(b, left_nodes)

    # Swap x and y coordinates to make vertical stacks horizontal
    for node in pos:
        pos[node] = (pos[node][1], -pos[node][0])  # Swap x and y coordinates

    # draw the graph
    plt.figure(figsize=(16, 8))
    nx.draw(b, pos, with_labels=True,
            node_color=['lightblue' if node in left_nodes else 'lightcoral' for node in b.nodes()],
            edge_color='gray', node_size=4000, font_size=9)
    plt.show()  # (5) done


from my_lib import explode_pathways


def draw_pathway_interactions_histogram(pathways: pd.DataFrame) -> None:
    boom = explode_pathways(pathways)

    # Plot
    plt.figure(figsize=(14, 6))

    # Create bin edges that align with the data
    bin_edges = np.arange(len(boom['drug'].unique()) + 1) - 0.5  # magic to somewhat align the x ticks

    # Plot histogram with aligned bins
    plt.hist(boom, bins=bin_edges)

    # Customize plot
    plt.xlabel('Drug Name')
    plt.ylabel('Count of Associated Pathways')
    plt.title('Number of Interaction Pathways per Drug')

    # Set x-ticks at the center of each bin
    plt.xticks(range(len(boom['drug'].unique())), boom['drug'].unique(), rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Show plot
    plt.show()  # (6) done


# code length is the same, imo easier code, and ticks are better aligned, but it's not a histogram :c
def draw_pathway_interactions_bar_chart(pathways: pd.DataFrame) -> None:
    boom = explode_pathways(pathways)

    # firstly we aggregate the data ourselves

    pathway_count = boom['drug'].value_counts()
    pathway_count.sort_index(inplace=True)

    # then we plot it as a bar chart
    plt.figure(figsize=(14, 6))

    plt.bar(pathway_count.index, pathway_count.values, color='skyblue')

    # Customize plot
    plt.xlabel('Drug Name')
    plt.ylabel('Count of Associated Pathways')
    plt.title('Number of Interaction Pathways per Drug')
    plt.xticks(rotation=45, ha='right')  # Rotate labels for readability -> align them to they still match
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Show plot
    plt.show()  # (6) done better ?


def draw_pie_chart(targets: pd.DataFrame) -> None:
    # repeat the last step so that if we run this cell more than once it doesn't mess itself up (useful for testing)
    locations = targets['cellular location'].value_counts()

    # Define a threshold for grouping small slices
    threshold = 1  # You can adjust this value (e.g., 1%, 2%, etc.)

    # Group small slices into "Other"
    small_slices = locations[locations / locations.sum() * 100 < threshold]
    locations = locations[locations / locations.sum() * 100 >= threshold]
    locations['Other'] = small_slices.sum()

    threshold = 20
    # Create an explosion array
    explode = [0.3 if (value / locations.sum() * 100) < threshold else 0 for value in locations]

    # Custom function to display percentages only if they are greater than 3%
    def show_greater_than_3(pct):
        return f'{pct:.1f}%' if pct > 2 else ''

    # nice colors
    colors = ['#001CF0', '#0038E2', '#0055D4', '#0071C6', '#008DB8', '#00AAAA',
              '#00C69C', '#00E28E', '#00FF80', ]

    # Plotting the pie chart
    plt.figure(figsize=(8, 8))  # Set the size of the chart
    plt.pie(
        locations,
        labels=locations.index,  # Use the unique locations as labels
        autopct=show_greater_than_3,
        # colors=plt.cm.Pastel1.colors,  # Use a pastel color scheme
        colors=colors,
        wedgeprops={'edgecolor': 'black', 'linewidth': 1.5},  # Add edges to the slices
        explode=explode
    )

    # Add a title
    plt.title('Distribution of Cellular Locations', fontsize=16, pad=20)

    # Equal aspect ratio ensures that the pie is drawn as a circle
    plt.axis('equal')

    # Display the chart
    plt.show()  # (8) done


def draw_summary_pie_chart(summary: pd.DataFrame) -> None:
    # Function to show actual values
    def absolute_value(val):
        a = round(val / 100 * sum(summary['number of drugs']))
        return f"{a}"  # Return as string

    # Plotting the pie chart
    plt.figure(figsize=(8, 8))  # Set the size of the chart
    plt.pie(
        summary['number of drugs'],
        labels=summary['status'],
        colors=plt.cm.Pastel1.colors,  # Use a pastel color scheme
        autopct=absolute_value,
    )

    # Add a title
    plt.title('Distribution of Cellular Locations', fontsize=16, pad=20)

    # Equal aspect ratio ensures that the pie is drawn as a circle
    plt.axis('equal')
    plt.legend(loc='upper right')

    # Display the chart
    plt.show()


def draw_gene_relations(gene_name: str, targets: pd.DataFrame, products: pd.DataFrame, drugs: pd.DataFrame) -> None:
    attackers = targets[targets['gene name'] == gene_name]['drug id'].copy()

    attacker_products = products[products['id'].isin(attackers)].copy()
    attacker_products.drop_duplicates(subset=['id', 'name'], inplace=True)
    attacker_products.reset_index(drop=True, inplace=True)
    attacker_products = attacker_products[['id', 'name']]
    attacker_products.loc[:, 'id'] = attacker_products['id'].map(drugs['name'])

    attackers_list = list(attackers.copy().map(drugs['name']))
    product_list = list(attacker_products['name'].copy())

    attackers_dict = {name: i + 1 for i, name in enumerate(attackers_list)}
    product_dict = {name: i + len(attackers_list) + 1 for i, name in enumerate(product_list)}

    attacker_products["id"] = attacker_products["id"].map(attackers_dict)
    attacker_products["name"] = attacker_products["name"].map(product_dict)

    def wrap_labels(label, width=12):
        return "\n".join(wrap(label, width, break_long_words=False))

    # Create a graph
    G = nx.Graph()

    node_number = 0

    # Add nodes
    G.add_node(node_number, label=gene_name, type="gene")

    for a in attackers_list:
        node_number += 1
        G.add_node(node_number, label=wrap_labels(a), type="drug")
    for p in product_list:
        node_number += 1
        G.add_node(node_number, label=wrap_labels(p), type="product")

    # Add edges with relationship labels
    edges = {
        "gene_target": [(0, i + 1) for i in range(len(attackers_list))],
        "contains_drug": list(attacker_products.itertuples(index=False, name=None))
    }

    # Define edge colors
    edge_colors = {
        "gene_target": "#66b3ff",  # "red",
        "contains_drug": "#99cc99"  # "green"
    }

    # Define node colors
    color_map = []
    for node in G:
        if G.nodes[node]["type"] == "gene":
            # color_map.append("red")
            color_map.append("#9370DB")  # purple
        elif G.nodes[node]["type"] == "drug":
            # color_map.append("blue")
            color_map.append("#66b3ff")  # light blue
        else:
            # color_map.append("green")
            color_map.append("#99cc99")  # light green

    # Add edges to the graph
    for key, edge_list in edges.items():
        G.add_edges_from(edge_list)

    pos = nx.spring_layout(G, seed=42)

    # Pobranie etykiet do wyświetlenia
    labels = nx.get_node_attributes(G, "label")

    # Draw nodes
    plt.figure(figsize=(16, 12))
    nx.draw(G, pos, with_labels=True, labels=labels, node_color=color_map, node_size=4000, edge_color='gray',
            font_size=10)

    # Draw edges with specific colors
    for relation, edge_list in edges.items():
        nx.draw_networkx_edges(G, pos, edgelist=edge_list, edge_color=edge_colors[relation], width=2)

    plt.show()


def draw_interactive_price_plot(prices: pd.DataFrame, x_grid:bool=False, y_grid:bool=False,
                                logarithmic:bool=True, scale:float=1.0) -> None:
    fig = px.scatter(prices, x="amount", y="cost",
                     color="unit",  # Kolorowanie według typu jednostki
                     title="Interactive Scatter Plot of Drug Prices",
                     labels={"amount": "Amount", "cost": "Price", "unit": "Unit Type"},
                     hover_name="description")

    if logarithmic:
        fig.update_yaxes(type="log", title="Price in USD (Log Scale)", showgrid=y_grid)
        fig.update_xaxes(type="log", title="Amount in mg (Log Scale)", showgrid=x_grid)
    else:
        fig.update_yaxes(title="Price in USD (Log Scale)", showgrid=y_grid)
        fig.update_xaxes(title="Amount in mg (Log Scale)", showgrid=x_grid)

    fig.update_layout(width=int(1600 * scale), height=int(1000 * scale))

    fig.update_traces(marker=dict(size=10, opacity=0.7))
    fig.show()
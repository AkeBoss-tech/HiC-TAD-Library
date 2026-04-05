import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def render_genome_scroller(genes_df: pd.DataFrame, current_pos: int, window_size: int = 1_000_000, deletion_region: tuple = None):
    """
    Render an interactive Plotly genome browser scroller with deletion overlay support.

    Args:
        genes_df: DataFrame with ['Feature', 'Start', 'End', 'Name', 'Type']
        current_pos: The center of the current view.
        window_size: Width of the scroller view.
        deletion_region: Tuple of (start, end) for a simulated deletion.
    """
    st.subheader("Interactive Genome Browser & Modification Planner")

    # Filter genes in view
    view_start = current_pos - (window_size // 2)
    view_end = current_pos + (window_size // 2)

    in_view = genes_df[
        (genes_df['End'] >= view_start) & (genes_df['Start'] <= view_end)
    ].copy()

    # Assign Y-levels to prevent overlaps
    if not in_view.empty:
        in_view = in_view.sort_values('Start')
        levels = []
        last_ends = []

        for _, row in in_view.iterrows():
            placed = False
            for i, last_end in enumerate(last_ends):
                if row['Start'] > last_end + (window_size * 0.05): # small gap
                    levels.append(i + 1)
                    last_ends[i] = row['End']
                    placed = True
                    break
            if not placed:
                levels.append(len(last_ends) + 1)
                last_ends.append(row['End'])

        in_view['Level'] = levels

    # Build Plotly Figure
    fig = go.Figure()

    # Add shadow for the deletion region
    if deletion_region:
        d_start, d_end = deletion_region
        fig.add_vrect(
            x0=d_start, x1=d_end,
            fillcolor="red", opacity=0.2,
            layer="below", line_width=0,
            annotation_text="Planned Deletion", 
            annotation_position="top left"
        )

    # Add features
    if not in_view.empty:
        for _, row in in_view.iterrows():
            color = 'blue' if row['Type'] == 'Gene' else 'orange' if row['Type'] == 'Enhancer' else 'green'
            
            fig.add_trace(go.Scatter(
                x=[row['Start'], row['End']],
                y=[row['Level'], row['Level']],
                mode='lines+markers',
                line=dict(color=color, width=10),
                name=row['Name'],
                text=f"{row['Name']} ({row['Type']})<br>{row['Start']}-{row['End']}",
                hoverinfo='text'
            ))

            # Label
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=row['Level'] + 0.3,
                text=row['Name'],
                showarrow=False,
                font=dict(size=10)
            )

    # Current focus indicator
    fig.add_vline(x=current_pos, line_width=2, line_dash="solid", line_color="black")

    fig.update_layout(
        height=450,
        showlegend=False,
        xaxis_title="Genomic Coordinates (hg38)",
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        xaxis=dict(range=[view_start, view_end], tickformat=",.0s"),
        margin=dict(l=20, r=20, t=20, b=20),
        hovermode='closest'
    )

    st.plotly_chart(fig, use_container_width=True)

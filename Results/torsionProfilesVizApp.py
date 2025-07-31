from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import pickle

def LoadData():
    with open('torsionHistogramsAllEnvs.pkl', 'rb') as f:
        data = pickle.load(f)
    return data

DATA = LoadData()

app = Dash()

# Requires Dash 2.17.0 or later
app.layout = [
    html.H1(children='Torsion Profiles', style={'font-family':'Arial', 'textAlign': 'center', 'marginBottom': '70px'}),
    html.Div(children='Entry number of the SMARTS pattern:', style={'font-family':'Arial', 'marginBottom': '10px'}),
    dcc.Input(id='entrynr', type='number', placeholder='entrynr', step=1, value=1),
    html.Div(children='SMARTS pattern:', style={'font-family':'Arial','marginBottom':'10px','marginTop':'10px'}),
    dcc.Textarea(
        id='smarts-info',
        style={'width': '100%', 'height': '25px'},
        readOnly=True
    ),
    dcc.Graph(id='graph')
]

@callback(
    Output('graph', 'figure'),
    [Input('entrynr', 'value')]
)
def update_graph(entrynr):
    if entrynr not in DATA:
        # Optionally, return an empty figure or a figure with an error message
        return go.Figure()
    # Create a 2x2 grid of subplots
    fig = make_subplots(rows=2, cols=2, subplot_titles=("Crystal", "Vacuum", "Water", "Hexane"))

    xHist = DATA[entrynr]['bins']
    binSize = 2 * np.pi / (36 * 2)
    crystalHist = DATA[entrynr]['crystal']
    vacHist = DATA[entrynr]['vac']
    tip3pHist = DATA[entrynr]['water']
    hexaneHist = DATA[entrynr]['hexane']

    # Find the maximum y value across all histograms
    ymax = max(np.max(crystalHist), np.max(vacHist), np.max(tip3pHist), np.max(hexaneHist))

    fig.update_layout(
        height=1000,
        width=1400,
        showlegend=True,
        xaxis_title="Angle (rad)",
        yaxis_title="Density",
        xaxis2_title="Angle (rad)",
        yaxis2_title="Density",
        xaxis3_title="Angle (rad)",
        yaxis3_title="Density",
        xaxis4_title="Angle (rad)",
        yaxis4_title="Density",
        barmode="overlay",
        legend=dict(title="Legend"),
        template="plotly_white"
    )

    fig.add_trace(go.Bar(
        x=xHist,
        y=crystalHist,
        width=binSize,
        marker_color="#C3DFDE",
        opacity=0.8,
        name="crystal"
    ), row=1, col=1)

    fig.add_trace(go.Bar(
        x=xHist,
        y=vacHist,
        width=binSize,
        marker_color="#DCB261",
        opacity=0.8,
        name="vac"
    ), row=1, col=2)

    fig.add_trace(go.Bar(
        x=xHist,
        y=tip3pHist,
        width=binSize,
        marker_color="#66A1BC",
        opacity=0.8,
        name="tip3p"
    ), row=2, col=1)

    fig.add_trace(go.Bar(
        x=xHist,
        y=hexaneHist,
        width=binSize,
        marker_color="#193D40",
        opacity=0.8,
        name="hexane"
    ), row=2, col=2)

    # Set the same ymax for all subplots
    fig.update_yaxes(range=[0, ymax * 1.05], row=1, col=1)
    fig.update_yaxes(range=[0, ymax * 1.05], row=1, col=2)
    fig.update_yaxes(range=[0, ymax * 1.05], row=2, col=1)
    fig.update_yaxes(range=[0, ymax * 1.05], row=2, col=2)

    # Update layout for subplots
    fig.update_layout(
        xaxis_title="Angle (rad)",
        yaxis_title="Density",
        barmode="overlay",
        legend=dict(title="Legend"),
        template="plotly_white"
    )
    return fig

@callback(
    Output('smarts-info', 'value'),
    [Input('entrynr', 'value')]
)
def update_smarts_info(entrynr):
    if entrynr not in DATA:
        return "No data available for this entry number."
    smarts = DATA[entrynr]['smarts']
    return f"{smarts}"

if __name__ == '__main__':
    app.run(debug=True)

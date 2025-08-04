import sys
sys.path.append('./Scripts')
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from Scripts import RetrieveRawDataGNN, RetrieveRawDataCsdModified, GetSmarts
from plotly.subplots import make_subplots

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
    html.Div(children='Energy threshold for the vac, tip3p and hexane environments:', style={'font-family':'Arial','marginTop': '20px', 'marginBottom': '10px'}),
    html.Div([dcc.Slider(
        id='threshold',
        min=5,
        max=50,
        step=5,
        value=25,  # Default value
        marks={i: str(i) for i in range(5, 55, 5)}
        )]),
        #,style= {'transform': 'scale(0.8)'}),
    dcc.Graph(id='graph')
]

@callback(
    Output('graph', 'figure'),
    [Input('entrynr', 'value'),
     Input('threshold', 'value')]
)
def update_graph(entrynr, threshold):
    # Create a 2x2 grid of subplots
    fig = make_subplots(rows=2, cols=2, subplot_titles=("Crystal", "Vacuum", "Water", "Hexane"))

    # retrieve data for all environments
    crystalData = RetrieveRawDataCsdModified(entrynr)
    vacData = RetrieveRawDataGNN(entrynr, threshold, 'vac')
    tip3pData = RetrieveRawDataGNN(entrynr, threshold, 'tip3p')
    hexaneData = RetrieveRawDataGNN(entrynr, threshold, 'hexane')

    # symmetrize data
    crystalData = np.concatenate((crystalData, 2*np.pi-crystalData))
    vacData = np.concatenate((vacData, 2*np.pi-vacData))
    tip3pData = np.concatenate((tip3pData, 2*np.pi-tip3pData))
    hexaneData = np.concatenate((hexaneData, 2*np.pi-hexaneData))

    binSize = 2 * np.pi / (36 * 2)
    bins = np.arange(0, 2 * np.pi + binSize, binSize)
    crystalHist, xHist = np.histogram(crystalData, bins=bins, density=True)
    vacHist, _ = np.histogram(vacData, bins=bins, density=True)
    tip3pHist, _ = np.histogram(tip3pData, bins=bins, density=True)
    hexaneHist, _ = np.histogram(hexaneData, bins=bins, density=True)

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
        x=xHist[:-1] + 0.5 * binSize,
        y=crystalHist,
        width=binSize,
        marker_color="#C3DFDE",
        opacity=0.8,
        name="crystal"
    ), row=1, col=1)

    fig.add_trace(go.Bar(
        x=xHist[:-1] + 0.5 * binSize,
        y=vacHist,
        width=binSize,
        marker_color="#DCB261",
        opacity=0.8,
        name="vac"
    ), row=1, col=2)

    fig.add_trace(go.Bar(
        x=xHist[:-1] + 0.5 * binSize,
        y=tip3pHist,
        width=binSize,
        marker_color="#66A1BC",
        opacity=0.8,
        name="tip3p"
    ), row=2, col=1)

    fig.add_trace(go.Bar(
        x=xHist[:-1] + 0.5 * binSize,
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
    smarts = GetSmarts(entrynr)
    return f"{smarts}"

if __name__ == '__main__':
    app.run(debug=True)

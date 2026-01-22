#!/usr/bin/env just --justfile

default: simulation mtb_tf

# Simulation Commands
simulation:
    pixi run simulation

simulation_create_model:
    pixi run create_model

simulation_find_metabolite_networks:
    pixi run find_metabolite_networks

simulation_metabolic_network_analysis:
    pixi run metabolic_network_analysis

simulation_ko_divergence:
    pixi run ko_divergence

simulation_mi_network:
    pixi run mi_network

simulation_density:
    pixi run density

simulation_imat:
    pixi run imat

simulation_viz:
    pixi run viz

# Mtb Transcription Factor Analysis
mtb_tf:
    pixi run mtb_tf

tf_reaction_info:
    pixi run tf_reaction_info

tf_metabolite_info:
    pixi run tf_metabolite_info

tf_gen_model:
    pixi run tf_gen_model

tf_sample_model:
    pixi run tf_sample_model

tf_divergence:
    pixi run tf_divergence

tf_ko_divergence:
    pixi run tf_ko_divergence

tf_target_density:
    pixi run tf_target_density

tf_reaction_centrality:
    pixi run tf_reaction_centrality

tf_metabolite_gsva:
    pixi run tf_metabolite_gsva

tf_viz:
    pixi run tf_viz

#!/usr/bin/env just --justfile

default: simulation mtb_tf

# Simulation Commands
simulation:
    pixi run simulation

create_model:
    pixi run create_model

find_metabolite_networks:
    pixi run find_metabolite_networks

metabolic_network_analysis:
    pixi run metabolic_network_analysis

ko_divergence:
    pixi run ko_divergence

mi_network:
    pixi run mi_network

density:
    pixi run density

imat:
    pixi run imat

viz:
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

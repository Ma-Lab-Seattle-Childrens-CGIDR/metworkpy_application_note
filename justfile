#!/usr/bin/env just --justfile

default: simulation

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

simulation:
    pixi run simulation

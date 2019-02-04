clc
clear all
close all

DATA = load('plasma-255.mat');
rgb = DATA.RGB;
min2(peaks)
hfig = figure()
hax = surf(peaks-min2(peaks))
view(-28,65)
colormap(hfig,rgb)
colorbar


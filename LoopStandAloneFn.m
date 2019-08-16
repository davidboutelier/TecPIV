clear all
close all
clc

R = [2:1:300]

for k =1:length(R)
    F = R(k)
    StandAlonePlotFn(F)
end
clear all; close all;

options = optimset('Display','iter');
x = fminbnd(@minimizer, 0, 100000,options)
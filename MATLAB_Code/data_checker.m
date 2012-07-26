clear all; close all; clc ;

load('oldDataMatrix.txt') ;
A = sparse(oldDataMatrix) ;
figure(1)
spy(A)

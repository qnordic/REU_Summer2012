clear all; close all; clc ;

load('oldDataMatrix.txt') ;
A = sparse(oldDataMatrix) ;
figure(1)
spy(A)

load('oldDataMatrixUpdate.txt') ;
B = sparse(oldDataMatrixUpdate) ;
figure(2)
spy(B)
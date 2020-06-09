clc
clear all
close all

Path = 'C:/Users/Rustam/Documents/Visual Studio 2017/Projects/ÑonvolutionalIntegralEquation/ÑonvolutionalIntegralEquation/Debug/';
ext ='.txt';

ns = load('C:/Users/Rustam/Documents/Visual Studio 2017/Projects/ÑonvolutionalIntegralEquation/ÑonvolutionalIntegralEquation/Debug/ParallelSeidelNS.txt') ;
ns_size = size(ns,2);
solverTypeParallel = 'Parallel';
solverTypeSequential = 'Sequential';

solveTypeSeidel = 'Seidel';
solveTypeReverse='Reverse';
CONSUMPTION = 'CONSUMPTION';
TIMES = 'TIMES';
ERROR ='ERROR';
NS = 'NS';

nsParSFile = strcat(Path, solverTypeParallel, solveTypeSeidel ,NS, ext);
errorParSFile = strcat(Path, solverTypeParallel, solveTypeSeidel ,ERROR, ext);
nsParS=load(nsParSFile);
errorParS = load(errorParSFile);
parallelSLabel = strcat(solverTypeParallel, solveTypeSeidel);
figure;
hold on
plot(nsParS,errorParS, 'r', 'Linewidth',1.3);
title('Error and N');
ylabel('Error');
xlabel('N');

nsParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,NS, ext);
errorParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,ERROR, ext);
nsParR=load(nsParRFile);
errorParR = load(errorParRFile);
parallelRLabel = strcat(solverTypeParallel, solveTypeReverse);
plot(nsParR,errorParR, 'k', 'Linewidth',1.3);
ylabel('Error');
xlabel('N');
    
hold on
nsSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,NS, ext);
errorSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,ERROR, ext);
nsSeqS=load(nsSeqSFile);
errorSeqS = load(errorSeqSFile);
seqSLabel = strcat(solverTypeSequential, solveTypeSeidel);
plot(nsSeqS,errorSeqS, 'g', 'Linewidth',1.3);
ylabel('Error');
xlabel('N');
    
hold on
nsSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,NS, ext);
errorSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,ERROR, ext);
nsSeqR=load(nsSeqRFile);
errorSeqR = load(errorSeqRFile);
seqRLabel = strcat(solverTypeSequential, solveTypeReverse);
plot(nsSeqR,errorSeqR, 'y', 'Linewidth',1.3);
ylabel('Error');
xlabel('N');

legend(parallelSLabel,parallelRLabel, seqSLabel, seqRLabel);

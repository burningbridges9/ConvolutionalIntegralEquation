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
timeParSFile = strcat(Path, solverTypeParallel, solveTypeSeidel ,TIMES, ext);
nsParS=load(nsParSFile);
timeParS = load(timeParSFile);
parallelSLabel = strcat(solverTypeParallel, solveTypeSeidel);
figure;
hold on
plot(nsParS,timeParS, 'r', 'Linewidth',1.3);
title('Time and N');
ylabel('Time');
xlabel('N');


nsParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,NS, ext);
timeParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,TIMES, ext);
nsParR=load(nsParRFile);
timeParR = load(timeParRFile);
parallelRLabel = strcat(solverTypeParallel, solveTypeReverse);
plot(nsParR,timeParR, 'k', 'Linewidth',1.3);
ylabel('Time');
xlabel('N');

    
hold on
nsSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,NS, ext);
timeSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,TIMES, ext);
nsSeqS=load(nsSeqSFile);
timeSeqS = load(timeSeqSFile);
seqSLabel = strcat(solverTypeSequential, solveTypeSeidel);
plot(nsSeqS,timeSeqS, 'g', 'Linewidth',1.3);
ylabel('Time');
xlabel('N');
    
hold on
nsSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,NS, ext);
timeSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,TIMES, ext);
nsSeqR=load(nsSeqRFile);
timeSeqR = load(timeSeqRFile);
seqRLabel = strcat(solverTypeSequential, solveTypeReverse);
plot(nsSeqR,timeSeqR, 'y', 'Linewidth',1.3);
ylabel('Time');
xlabel('N');

legend(parallelSLabel,parallelRLabel, seqSLabel, seqRLabel);

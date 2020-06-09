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

for i = 1 : ns_size
    N = int2str(ns(i));
    consParSFile = strcat(Path, solverTypeParallel, solveTypeSeidel ,N, CONSUMPTION, ext);
    timesParSFile = strcat(Path, solverTypeParallel, solveTypeSeidel ,N, TIMES, ext);
    consParS=load(consParSFile) * 24.0 * 3600.0;
    timesParS = load(timesParSFile) / 3600.0;
    parallelSLabel = strcat('N = ', N, ' ',solverTypeParallel,' ', solveTypeSeidel);
    figure;
    plot(timesParS,consParS, 'r', 'Linewidth',1.2);
    labelForSolve = strcat('Solve at N = ', N);
    title(labelForSolve);
    ylabel('Consumption');
    xlabel('Time');
    hold on
    
    consParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,N, CONSUMPTION, ext);
    timesParRFile = strcat(Path, solverTypeParallel, solveTypeReverse ,N, TIMES, ext);
    consParR=load(consParRFile) * 24.0 * 3600.0;
    timesParR = load(timesParRFile) / 3600.0;
    parallelRLabel = strcat('N = ', N, ' ',solverTypeParallel,' ', solveTypeReverse);
    plot(timesParR,consParR, 'k', 'Linewidth',1.2);
    ylabel('Consumption');
    xlabel('Time');
    hold on
    
    consSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,N, CONSUMPTION, ext);
    timesSeqSFile = strcat(Path, solverTypeSequential, solveTypeSeidel ,N, TIMES, ext);
    consSeqS=load(consSeqSFile) * 24.0 * 3600.0;
    timesSeqS = load(timesSeqSFile) / 3600.0;
    seqSLabel = strcat('N = ', N, ' ',solverTypeSequential,' ', solveTypeSeidel);
    plot(timesSeqS,consSeqS, 'g', 'Linewidth',1.2);
    ylabel('Consumption');
    xlabel('Time');
    hold on
    
    consSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,N, CONSUMPTION, ext);
    timesSeqRFile = strcat(Path, solverTypeSequential, solveTypeReverse ,N, TIMES, ext);
    consSeqR=load(consSeqRFile) * 24.0 * 3600.0;
    timesSeqR = load(timesSeqRFile) / 3600.0;
    seqRLabel = strcat('N = ', N, ' ',solverTypeSequential,' ', solveTypeReverse);
    plot(timesSeqR,consSeqR, 'y', 'Linewidth',1.2);
    ylabel('Consumption');
    xlabel('Time');
    hold on
    
    legend(parallelSLabel,parallelRLabel,seqSLabel, seqRLabel);
    hold off
end

%% Model Predictive Control
% MPC for ACC System
clc;
clear all;
close all;

%%display
disp('--- Model Predictive Control of a Adaptive Cruise Control ---');
disp('Choose one of the options below');
disp('1. ACC');
disp('2. ACC MPC and ACC LQR Comparison');
disp('3. MPC of ACC with varying prediction horizon N');
disp('4. MPC of ACC with varying Q');
disp('5. MPC of ACC with varying R');
disp('0. Exit');
pause(1);
choice = input('Please choose --> ');
switch choice
    case 1
        close all;
        figure(1);
        ACC_MPC_Final;
        disp('Do you want to run another simulation?');
        opt = input('Enter your option [Y/any other key] --> ','s');
        if opt == 'Y'
            main;
        else
            disp('Goodbye!');
        end
    case 2
        close all;
        figure(2);
        LQR_ACC;
        disp('Do you want to run another simulation?');
        opt = input('Enter your option [Y/any other key] --> ','s');
        if opt == 'Y'
            main;
        else
            disp('Goodbye!');
        end
    case 3
        close all;
        figure(3);
        ACC_MPC_varying_N;
        disp('Do you want to run another simulation?');
        opt = input('Enter your option [Y/any other key] --> ','s');
        if opt == 'Y'
            main;
        else
            disp('Goodbye!');
        end
    case 4
        close all;
        figure(4);
        ACC_MPC_varying_Q;
        disp('Do you want to run another simulation?');
        opt = input('Enter your option [Y/any other key] --> ','s');
        if opt == 'Y'
            main;
        else
            disp('Goodbye!');
        end
    case 5
        close all;
        figure(5);
        ACC_MPC_varying_R;
        disp('Do you want to run another simulation?');
        opt = input('Enter your option [Y/any other key] --> ','s');
        if opt == 'Y'
            main;
        else
            disp('Goodbye!');
        end
    case 0
        disp('Goodbye!');
    otherwise 
        disp('Please input a valid number.');
        pause(2);
        main;
end

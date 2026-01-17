% RUN_ME - Quick Start Guide for Power Quality Analysis Project
%
% This script provides a menu-driven interface to run all analysis options
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025

clear all; close all; clc;

fprintf('\n');
fprintf('================================================================================\n');
fprintf('     POWER QUALITY ANALYSIS USING DTFS - PROJECT MENU\n');
fprintf('================================================================================\n');
fprintf('Author: Manikanta Gonugondla (230906450)\n');
fprintf('Course: Digital Signal Processing (FISAC Assessment)\n');
fprintf('Institution: MIT & IIT Madras\n');
fprintf('================================================================================\n\n');

fprintf('Available Analysis Options:\n\n');

fprintf('COMPREHENSIVE ANALYSIS:\n');
fprintf('  [1] Complete Analysis with ALL Visualizations (7 figures)\n');
fprintf('      - Time domain, frequency domain, phase spectrum\n');
fprintf('      - Harmonic pie charts, THD analysis\n');
fprintf('      - Filtering results, power quality dashboard\n');
fprintf('      File: comprehensive_analysis.m\n\n');

fprintf('  [2] Main Power Quality Analysis (4 figures)\n');
fprintf('      - Original main analysis script\n');
fprintf('      - Time/frequency comparison, filtering, metrics\n');
fprintf('      File: main_power_quality_analysis.m\n\n');

fprintf('INDIVIDUAL SCENARIO ANALYSIS:\n');
fprintf('  [3] LED Lighting Scenario\n');
fprintf('      - 3rd harmonic dominant analysis\n');
fprintf('      File: scenarios/scenario_led_lighting.m\n\n');

fprintf('  [4] Motor Drive (VFD) Scenario\n');
fprintf('      - 5th & 7th harmonic analysis with filtering\n');
fprintf('      File: scenarios/scenario_motor_drive.m\n\n');

fprintf('  [5] Data Center (SMPS) Scenario\n');
fprintf('      - Multiple odd harmonics analysis\n');
fprintf('      File: scenarios/scenario_data_center.m\n\n');

fprintf('VALIDATION TESTS:\n');
fprintf('  [6] Run All Test Scripts\n');
fprintf('      - Validates all 6 core functions\n\n');

fprintf('  [7] Individual Test Scripts\n');
fprintf('      - test_calculate_dtfs.m\n');
fprintf('      - test_synthesize_dtfs.m\n');
fprintf('      - test_generate_signals.m\n');
fprintf('      - test_calculate_thd.m\n');
fprintf('      - test_power_quality_metrics.m\n');
fprintf('      - test_harmonic_filter.m\n\n');

fprintf('DOCUMENTATION:\n');
fprintf('  [8] Compile LaTeX Report (requires pdflatex)\n');
fprintf('      File: report.tex\n\n');

fprintf('  [9] View Project Summary\n');
fprintf('      File: PROJECT_SUMMARY.txt\n\n');

fprintf('  [0] Exit\n\n');

fprintf('================================================================================\n');
choice = input('Enter your choice [1-9, 0 to exit]: ', 's');
fprintf('================================================================================\n\n');

switch choice
    case '1'
        fprintf('Running COMPREHENSIVE ANALYSIS...\n');
        fprintf('This will generate 7 comprehensive figures.\n');
        fprintf('Press Ctrl+C to cancel, or wait 3 seconds to continue...\n');
        pause(3);
        comprehensive_analysis;

    case '2'
        fprintf('Running MAIN POWER QUALITY ANALYSIS...\n');
        fprintf('This will generate 4 figures.\n');
        fprintf('Press Ctrl+C to cancel, or wait 3 seconds to continue...\n');
        pause(3);
        main_power_quality_analysis;

    case '3'
        fprintf('Running LED LIGHTING SCENARIO...\n');
        cd scenarios;
        scenario_led_lighting;
        cd ..;

    case '4'
        fprintf('Running MOTOR DRIVE SCENARIO...\n');
        cd scenarios;
        scenario_motor_drive;
        cd ..;

    case '5'
        fprintf('Running DATA CENTER SCENARIO...\n');
        cd scenarios;
        scenario_data_center;
        cd ..;

    case '6'
        fprintf('Running ALL TEST SCRIPTS...\n\n');
        fprintf('Test 1/6: DTFS Calculation\n');
        test_calculate_dtfs;
        fprintf('\nTest 2/6: DTFS Synthesis\n');
        test_synthesize_dtfs;
        fprintf('\nTest 3/6: Signal Generation\n');
        test_generate_signals;
        fprintf('\nTest 4/6: THD Calculation\n');
        test_calculate_thd;
        fprintf('\nTest 5/6: Power Quality Metrics\n');
        test_power_quality_metrics;
        fprintf('\nTest 6/6: Harmonic Filter\n');
        test_harmonic_filter;
        fprintf('\nAll tests completed!\n');

    case '7'
        fprintf('Individual test scripts:\n');
        fprintf('  Run: test_calculate_dtfs\n');
        fprintf('  Run: test_synthesize_dtfs\n');
        fprintf('  Run: test_generate_signals\n');
        fprintf('  Run: test_calculate_thd\n');
        fprintf('  Run: test_power_quality_metrics\n');
        fprintf('  Run: test_harmonic_filter\n');

    case '8'
        fprintf('Compiling LaTeX report...\n');
        if ispc
            fprintf('On Windows, run: pdflatex report.tex\n');
            fprintf('Run it twice to resolve references.\n');
        elseif ismac || isunix
            fprintf('On Mac/Linux, run: ./compile_report.sh\n');
            fprintf('Or: pdflatex report.tex (run twice)\n');
            system('./compile_report.sh');
        end

    case '9'
        fprintf('Opening PROJECT_SUMMARY.txt...\n');
        if ispc
            system('type PROJECT_SUMMARY.txt');
        else
            system('cat PROJECT_SUMMARY.txt');
        end

    case '0'
        fprintf('Exiting...\n');
        return;

    otherwise
        fprintf('Invalid choice! Please run RUN_ME again and select 1-9 or 0.\n');
end

fprintf('\n');
fprintf('================================================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('================================================================================\n');
fprintf('\nTo run another analysis, execute: RUN_ME\n');
fprintf('To see all options, type: help RUN_ME\n\n');

fprintf('RECOMMENDED SEQUENCE FOR FIRST-TIME USERS:\n');
fprintf('  1. Run option [6] to validate all functions (tests)\n');
fprintf('  2. Run option [1] for comprehensive analysis (7 figures)\n');
fprintf('  3. Run individual scenarios [3-5] for detailed analysis\n');
fprintf('  4. Compile LaTeX report [8] for documentation\n\n');

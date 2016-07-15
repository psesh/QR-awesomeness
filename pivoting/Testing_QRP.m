% Comparing my version of MGS & Householder with pivoting
% with MATLAB's qr for random matrices
clear; close all; clc;
A = rand(90,60); [m,n] = size(A); % random matrix
[Q, R, P] = qr_MGS_pivoting(A);
[Q2, R2, P2] = qr_Householder_pivoting(A);
[Q3, R3, P3] = qr(A,'vector');

%%
figure1 = figure;
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; grid on;
plot(P, 'bo', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MGS'); 
plot(P2, 'ks', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'House'); 
plot(P3, 'rx', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MATLAB-qr'); 
legend show;
xlabel('Pivot columns'); hold off;

%%
figure2 = figure;
subplot(2,3,1);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(Q' * Q))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin')

subplot(2,3,2)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(Q2' * Q2))'); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 m]); ylim([1 m]);
view([90 90]);
title('Householder')

subplot(2,3,3)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(Q3' * Q3))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 m]); ylim([1 m]);
view([90 90]);
title('MATLAB')

subplot(2,3,4);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(R'))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin')

subplot(2,3,5);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(R2'))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 m]); ylim([1 n]); 
view([90 90]);
title('Householder')

subplot(2,3,6);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(R3'))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 m]); ylim([1 n]); 
view([90 90]);
title('MATLAB')

hold off;

%%
figure3 = figure;
subplot(1,3,1);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; 
pcolor(log10(abs(A * permutation(P') - Q * R))); shading flat;
caxis([-16 0]);  axis equal
xlabel('n'); ylabel('m');
xlim([1 n]); ylim([1 m]); 
view([90 90]);
title('MGS-Thin')

subplot(1,3,2)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on;
pcolor(log10(abs(A * permutation(P2') - Q2 * R2))); shading flat;
caxis([-16 0]);axis equal
xlabel('n'); ylabel('m');
xlim([1 n]); ylim([1 m]);
view([90 90]);
title('Householder')

subplot(1,3,3)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on;
pcolor(log10(abs(A * permutation(P3') - Q3 * R3))); shading flat;
caxis([-16 0]); axis equal
xlabel('n'); ylabel('m');
xlim([1 n]); ylim([1 m]);
view([90 90]);
title('MATLAB')

axes1 = axes('Visible','off','Parent',figure3,'Clipping','off',...
    'Position',[0 0 1 1]);
text('Parent',axes1,'VerticalAlignment','top',...
    'HorizontalAlignment','center',...
    'String','\bf (AP - QR)',...
    'Position',[0.525 0.904761904761905 0], 'FontSize', 18);
hold off;
% Comparing my version of MGS & Householder with pivoting
% with MATLAB's qr for random matrices
clear; close all; clc;
A = rand(90,60); [m,n] = size(A); % random matrix
[Q, R, P] = qr_MGS_pivoting(A);
[Q2, R2, P2] = qr_MGS_pivoting_custom(A, [1 2 3 4 5 6 7 8 9 10 11 12 13]);

% Compare pivots!
figure1 = figure;
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; grid on;
plot(P, 'bo', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MGS'); 
plot(P2, 'ks', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MGS-fixed'); 
legend show;
xlabel('Pivot columns'); hold off;

%%
figure2 = figure;
subplot(2,2,1);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(Q' * Q))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin')

subplot(2,2,2)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(Q2' * Q2))'); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin-Fixed')

subplot(2,2,3);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(R'))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin')

subplot(2,2,4);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; axis equal;
pcolor(log10(abs(R2'))); shading flat;
caxis([-5 0]);
xlabel('m'); ylabel('n');
xlim([1 min(n,m)]); ylim([1 min(n,m)]); 
view([90 90]);
title('MGS-Thin-Fixed')


hold off;

%%
figure3 = figure;
subplot(1,2,1);
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; 
pcolor(log10(abs(A * permutation(P') - Q * R))); shading flat;
caxis([-16 0]);  axis equal
xlabel('n'); ylabel('m');
xlim([1 n]); ylim([1 m]); 
view([90 90]);
title('MGS-Thin')

subplot(1,2,2)
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on;
pcolor(log10(abs(A * permutation(P2') - Q2 * R2))); shading flat;
caxis([-16 0]);axis equal
xlabel('n'); ylabel('m');
xlim([1 n]); ylim([1 m]); 
view([90 90]);
title('MGS-Thin-Fixed')


axes1 = axes('Visible','off','Parent',figure3,'Clipping','off',...
    'Position',[0 0 1 1]);
text('Parent',axes1,'VerticalAlignment','top',...
    'HorizontalAlignment','center',...
    'String','\bf (AP - QR)',...
    'Position',[0.525 0.904761904761905 0], 'FontSize', 18);
hold off;
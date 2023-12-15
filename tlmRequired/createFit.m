function [fitresult, gof] = createFit(xx, yy, a)
%CREATEFIT(XX,YY,A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xx
%      Y Input : yy
%      Z Output: a
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 13-May-2020 21:14:59 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( xx, yy, a );

% Set up fittype and options.
ft = fittype( 'poly11' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'a vs. xx, yy', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'xx', 'Interpreter', 'none' );
% ylabel( 'yy', 'Interpreter', 'none' );
% zlabel( 'a', 'Interpreter', 'none' );
% grid on
% view( -15.2, 4.9 );



%Inclass 22
%GB comments
1 100
2 100
overall 100

%1. Consider the case of the auto-activating gene that we discussed in class
%today. Make a bifurcation diagram for this system by varying the
%activated transcription rate for three cases - in which 1, 4, or 8 copies of the
%transcripton factor are necessary to activate transciption. Comment on your
%results. 

% Case 1 - 1 copy necessary
figure
hold on
ku = 0.1;
for kb = 0:0.05:5
    
    % Define polynomial of transcribed protein where rate of generation = 0
    coeffs = [1 (1-kb) -ku];
    rts = roots(coeffs);
    % Take only real roots
    rts_real = rts(imag(rts) == 0);
    plot(kb*ones(length(rts_real),1),rts_real,'r.')
    
end

% Label plot
xlabel('k_b')
ylabel('Fixed points')
title('n = 1')
set(gca, 'FontSize', 24)

% Case 2 - 4 copies necessary
figure
hold on
ku = 0.1;
for kb = 0:0.05:5
    
    % Define polynomial of transcribed protein where rate of generation = 0
    coeffs = [1 -kb 0 0 0 0 0 0 1 -ku];
    rts = roots(coeffs);
    % Take only real roots
    rts_real = rts(imag(rts) == 0);
    plot(kb*ones(length(rts_real),1),rts_real,'r.')
    
end

% Label plot
xlabel('k_b')
ylabel('Fixed points')
title('n = 4')
set(gca, 'FontSize', 24)

% Case 3 - 8 copies necessary
figure
hold on
ku = 0.1;
for kb = 0:0.05:5
    
    % Define polynomial of transcribed protein where rate of generation = 0
    coeffs = [1 -kb 1 -ku];
    rts = roots(coeffs);
    % Take only real roots
    rts_real = rts(imag(rts) == 0);
    plot(kb*ones(length(rts_real),1),rts_real,'r.')
    
end

% Label plot
xlabel('k_b')
ylabel('Fixed points')
title('n = 8')
set(gca, 'FontSize', 24)

% In case 1, where only 1 copy of a transcription factor is needed to
% initiate auto-activation, two stable solutions are present. One of these
% is non-physical, with a negative concentration of the transription
% factor. This appears because the rate of decay is artificially positive
% and balances the rate of production, which is artificially negative. The
% second stable solution is real and represents the case where production
% outweighs degradation (or dilution). This occurs at all possible kb
% values if x is large enough, since basal expression occurs.

% In case 2, where 4 copies of the transcription factor are needed to
% initiate activation, two stable solutions and one unstable solution are
% present. The case where the basal expression balances the degradation is
% the first stable solution (no effective activation). The second stable
% solution occurs where kb becomes large enough so that auto-activation and
% basal expression added together balance degradation. kb has become large
% enough to sustain auto-activation at higher x concentrations. The
% unstable solution is the case where auto-activation is occurring, but
% degradation has become large, so the rates of low auto-activation plus
% basal expression balance the rate of degradation. A perturbation of
% higher x will promote more auto-activation, driving the solution up
% towards the upper stable solution. A perturbation of lower x will cause
% degradation to dominate, pushing the solution down towards the lower
% stable solution with no auto-regulation.

% In case 3, where 8 copies of the transription factor are needed to
% initiate activation, the same trends are seen as in case 3. However, the
% unstable solution only exists in a much shorter range of kb values, and
% the lower stable solution disappears along with the unstable solution. At
% high kb values, auto-activation becomes so significant that it pushes all
% x concentrations up to the upper steady state. The rate of promotion
% becomes so high that it cannot be balanced by degradation.

% 2. Make a similar diagram for the case of an autorepressing gene in the
% case that 2 copies are need to turn off the gene. 

% Case 4 - 2 copies repress
figure
hold on
ku = 10;
for kb = 0:0.05:5
    
    % Define polynomial of transcribed protein where rate of generation = 0
    coeffs = [1 -kb 1 -ku];
    rts = roots(coeffs);
    % Take only real roots
    rts_real = rts(imag(rts) == 0);
    plot(kb*ones(length(rts_real),1),rts_real,'r.')
    
end

% Label plot
xlabel('k_b')
ylabel('Fixed points')
title('n = 2, repressed')
set(gca, 'FontSize', 24)

% In this case, the transcription factor is auto-repressed. One stable
% solution is present. At lower values of kb, auto-repression is more
% pronounced. x will auto-repress itself to a greater degree, as its rate
% of production in the bound state is much lower, so the stable
% concentration of x is relatively low. When kb becomes larger (closer to
% ku), auto-repression becomes less significant. Basal expression and
% repressed bound expression generates a higher concentration of x, because
% the repression does not lower the rate of production as much.


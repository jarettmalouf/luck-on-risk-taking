% Columns in 'data' matrices (other columns are blank):
% 1. Trial # for given participant
% 3. Certain points
% 4. Gain points
% 5. Loss points
% 7. Chosegamble (1/0)
% 8. Outcome (= Certain points, Gain points, or Loss points)
% 9. Response time

load('lucky_14122022.mat');
lucky = alldata;
load('unlucky_14122022.mat');
unlucky = alldata;

TRIAL_NO = 1;
CERTAIN = 3;
GAIN = 4;
LOSS = 5;
CHOSE_GAMBLE = 7;
OUTCOME = 8;
RT = 9;
RESULT = 15; %1 safe, 2 loss, 3 won

SAFE_RESULT = 1;
LOSS_RESULT = 2;
WON_RESULT = 3;

N_LUCKY = length(lucky);
N_UNLUCKY = length(unlucky);
T = length(lucky(1).data);
INIT_SCORE = 1000;
Data cleanup
% Overall gambling

for i = 1 : N_LUCKY
    individual_data = lucky(i).data;

    gamble_trials = individual_data(individual_data(:,CHOSE_GAMBLE) == 1);
    lucky(i).percent_gambles_taken = length(gamble_trials) / T;

    % Biased gambling (did they catch on?)
    gamble_trials = individual_data(individual_data(T/2+1:T,CHOSE_GAMBLE) == 1) + T/2;
    lucky(i).percent_gambles_taken_biased = length(gamble_trials) / (T/2);
end

for i = 1 : N_UNLUCKY
    individual_data = unlucky(i).data;

    gamble_trials = individual_data(individual_data(:,CHOSE_GAMBLE) == 1);
    unlucky(i).percent_gambles_taken = length(gamble_trials) / T;

    gamble_trials = individual_data(individual_data(T/2+1:T,CHOSE_GAMBLE) == 1) + T/2;
    unlucky(i).percent_gambles_taken_biased = length(gamble_trials) / (T/2);
end

GAMBLING_THRESHOLD = 0.2;
lucky = lucky([lucky.percent_gambles_taken] > GAMBLING_THRESHOLD & [lucky.percent_gambles_taken] < 1 - GAMBLING_THRESHOLD);
unlucky = unlucky([unlucky.percent_gambles_taken] > GAMBLING_THRESHOLD & [unlucky.percent_gambles_taken] < 1 - GAMBLING_THRESHOLD);

GAMBLING_THRESHOLD = 0.1;
lucky = lucky([lucky.percent_gambles_taken_biased] > GAMBLING_THRESHOLD & [lucky.percent_gambles_taken_biased] < 1 - GAMBLING_THRESHOLD);
unlucky = unlucky([unlucky.percent_gambles_taken_biased] > GAMBLING_THRESHOLD & [unlucky.percent_gambles_taken_biased] < 1 - GAMBLING_THRESHOLD);

N_LUCKY = length(lucky);
N_UNLUCKY = length(unlucky);

data = [lucky; unlucky];
N = length(data);
Verification of methodology
Let's start by verifying that the proportion of gambles in the first half of each trial were 50/50 and that they changed to 75/25 and 25/75 respectively in the lucky and unlucky groups.
for i = 1 : N_LUCKY

    individual_data = lucky(i).data;

    %     first half
    gamble_trials = individual_data(individual_data(1:T/2,CHOSE_GAMBLE) == 1);
    lucky(i).percent_won_gambles_fair = sum(individual_data(gamble_trials,RESULT) == WON_RESULT) / length(gamble_trials);

    %     second half
    gamble_trials = individual_data(individual_data(T/2 + 1:T,CHOSE_GAMBLE) == 1) + T/2;
    lucky(i).percent_won_gambles_biased = sum(individual_data(gamble_trials,RESULT) == WON_RESULT) / length(gamble_trials);

end

mean_percent_fair_won_lucky = mean([lucky.percent_won_gambles_fair]);
mean_percent_biased_won_lucky = mean([lucky.percent_won_gambles_biased]);

sem_percent_fair_won_lucky = sem([lucky.percent_won_gambles_fair]);
sem_percent_biased_won_lucky = sem([lucky.percent_won_gambles_biased]);

for i = 1 : N_UNLUCKY

    individual_data = unlucky(i).data;

    %     first half
    gamble_trials = individual_data(individual_data(1:T/2,CHOSE_GAMBLE) == 1);
    unlucky(i).percent_won_gambles_fair = sum(individual_data(gamble_trials,RESULT) == WON_RESULT) / length(gamble_trials);

    %     second half
    gamble_trials = individual_data(individual_data(T/2 + 1:T,CHOSE_GAMBLE) == 1) + T/2;
    unlucky(i).percent_won_gambles_biased = sum(individual_data(gamble_trials,RESULT) == WON_RESULT) / length(gamble_trials);
    
end

mean_percent_fair_won_unlucky = mean([unlucky.percent_won_gambles_fair]);
mean_percent_biased_won_unlucky = mean([unlucky.percent_won_gambles_biased]);

sem_percent_fair_won_unlucky = sem([unlucky.percent_won_gambles_fair]);
sem_percent_biased_won_unlucky = sem([unlucky.percent_won_gambles_biased]);

% Visualization

fair_means = [mean_percent_fair_won_lucky,  mean_percent_fair_won_unlucky];
biased_means = [mean_percent_biased_won_lucky,  mean_percent_biased_won_unlucky];

fair_sems = [sem_percent_fair_won_lucky,  sem_percent_fair_won_unlucky];
biased_sems = [sem_percent_biased_won_lucky,  sem_percent_biased_won_unlucky];

p_lucky_fair = signrank([lucky.percent_won_gambles_fair], .5);
p_lucky_biased = signrank([lucky.percent_won_gambles_biased], .75);
p_unlucky_fair = signrank([unlucky.percent_won_gambles_fair], .5);
p_unlucky_biased = signrank([unlucky.percent_won_gambles_biased], .25);

figure
labels = {'Lucky','Unlucky'};
b = bar([fair_means; biased_means]', 'grouped');
hold on; yline(.5,"b:")
hold on
set(gca,'xtick',1:length(labels),'xticklabel',labels);
x = [b(1).XEndPoints; b(2).XEndPoints];
errorbar([b(1).XEndPoints; b(2).XEndPoints]', [fair_means; biased_means]', [fair_sems; biased_sems]', '.');
title(sprintf("Online data (n = %d)", N))
ylim([0 1])
ylabel('Percent gambles won');
legend({'Fair', 'Biased'})
Distribution of luck
LUCK_THRESHOLD = .05;
bin_frequencies = zeros(length(0:LUCK_THRESHOLD:1),1);

for i = 1 : N

    individual_data = data(i).data;

    gamble_trials = individual_data(individual_data(:,CHOSE_GAMBLE) == 1);
    data(i).percent_won_gambles = sum(individual_data(gamble_trials,RESULT) == WON_RESULT) / length(gamble_trials);
    data(i).binned_percent_won_gambles = round(data(i).percent_won_gambles / LUCK_THRESHOLD) * LUCK_THRESHOLD;

    bin_index = round(1 + (data(i).binned_percent_won_gambles / LUCK_THRESHOLD));
    bin_frequencies(bin_index) = bin_frequencies(bin_index) + 1;

end

figure
bar(0:LUCK_THRESHOLD:1, bin_frequencies)
xlabel("Percent gambles won")
ylabel("Participant distribution")
xlim([0 1])
ylim([0 5])
title(sprintf("All participants (n = %d)", N))
Risk-taking across gamble domains
Next we'll want to investigate whether changing the probability of a successful risky outcome (i.e. the perception of luck) affects the participants' attitude toward risk. That is, within each group (lucky/unlucky) and within each domain (gain/mixed/loss), whether they were more or less willing to take risks relative to the unbiased first half of the experiment.
% LUCKY

for i = 1 : N_LUCKY

    individual_data = lucky(i).data;

    %     first half
    gain_trials = individual_data(individual_data(1:T/2,4) > 0 & individual_data(1:T/2,5) == 0);
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    lucky(i).percent_risk_in_gain_fair = num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(1:T/2,4) == 0 & individual_data(1:T/2,5) < 0);
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    lucky(i).percent_risk_in_loss_fair = num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(1:T/2,4) > 0 & individual_data(1:T/2,5) < 0);
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    lucky(i).percent_risk_in_mixed_fair = num_mixed_gambles / length(mixed_trials);

    %     second half
    gain_trials = individual_data(individual_data(T/2 + 1:T,4) > 0 & individual_data(T/2 + 1:T,5) == 0) + T/2;
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    lucky(i).percent_risk_in_gain_biased = num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(T/2 + 1:T,4) == 0 & individual_data(T/2 + 1:T,5) < 0) + T/2;
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    lucky(i).percent_risk_in_loss_biased = num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(T/2 + 1:T,4) > 0 & individual_data(T/2 + 1:T,5) < 0) + T/2;
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    lucky(i).percent_risk_in_mixed_biased = num_mixed_gambles / length(mixed_trials);

end

% Visualization

mean_percent_risk_in_gain_fair = nanmean([lucky.percent_risk_in_gain_fair]);
mean_percent_risk_in_loss_fair = nanmean([lucky.percent_risk_in_loss_fair]);
mean_percent_risk_in_mixed_fair = nanmean([lucky.percent_risk_in_mixed_fair]);

mean_percent_risk_in_gain_biased = nanmean([lucky.percent_risk_in_gain_biased]);
mean_percent_risk_in_loss_biased = nanmean([lucky.percent_risk_in_loss_biased]);
mean_percent_risk_in_mixed_biased = nanmean([lucky.percent_risk_in_mixed_biased]);

sem_percent_risk_in_gain_fair = sem([lucky.percent_risk_in_gain_fair]);
sem_percent_risk_in_loss_fair = sem([lucky.percent_risk_in_loss_fair]);
sem_percent_risk_in_mixed_fair = sem([lucky.percent_risk_in_mixed_fair]);

sem_percent_risk_in_gain_biased = sem([lucky.percent_risk_in_gain_biased]);
sem_percent_risk_in_loss_biased = sem([lucky.percent_risk_in_loss_biased]);
sem_percent_risk_in_mixed_biased = sem([lucky.percent_risk_in_mixed_biased]);

fair_percent_risk_means = [mean_percent_risk_in_gain_fair, mean_percent_risk_in_loss_fair, mean_percent_risk_in_mixed_fair]
biased_percent_risk_means = [mean_percent_risk_in_gain_biased,  mean_percent_risk_in_loss_biased, mean_percent_risk_in_mixed_biased]

fair_percent_risk_sems = [sem_percent_risk_in_gain_fair, sem_percent_risk_in_loss_fair, sem_percent_risk_in_mixed_fair];
biased_percent_risk_sems = [sem_percent_risk_in_gain_biased,  sem_percent_risk_in_loss_biased, sem_percent_risk_in_mixed_biased];

figure
subplot(1,2,1)
labels = {'Gain','Loss','Mixed'};
b = bar([fair_percent_risk_means; biased_percent_risk_means]', 'grouped');
% hold on; yline(.5,"b:")
hold on
set(gca,'xtick',1:length(labels),'xticklabel',labels);
errorbar([b(1).XEndPoints; b(2).XEndPoints]', [fair_percent_risk_means; biased_percent_risk_means]', ...
        [fair_percent_risk_sems; biased_percent_risk_sems]', '.');
title(sprintf("Lucky participants (n = %d)", N_LUCKY))
ylim([0 1])
ylabel('Percent risks taken');
legend({'Fair', 'Biased'})

% UNLUCKY

for i = 1 : N_UNLUCKY

    individual_data = unlucky(i).data;

    %     first half
    gain_trials = individual_data(individual_data(1:T/2,4) > 0 & individual_data(1:T/2,5) == 0);
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    unlucky(i).percent_risk_in_gain_fair = num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(1:T/2,4) == 0 & individual_data(1:T/2,5) < 0);
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    unlucky(i).percent_risk_in_loss_fair = num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(1:T/2,4) > 0 & individual_data(1:T/2,5) < 0);
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    unlucky(i).percent_risk_in_mixed_fair = num_mixed_gambles / length(mixed_trials);

    %     second half
    gain_trials = individual_data(individual_data(T/2 + 1:T,4) > 0 & individual_data(T/2 + 1:T,5) == 0) + T/2;
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    unlucky(i).percent_risk_in_gain_biased = num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(T/2 + 1:T,4) == 0 & individual_data(T/2 + 1:T,5) < 0) + T/2;
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    unlucky(i).percent_risk_in_loss_biased = num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(T/2 + 1:T,4) > 0 & individual_data(T/2 + 1:T,5) < 0) + T/2;
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    unlucky(i).percent_risk_in_mixed_biased = num_mixed_gambles / length(mixed_trials);

end

% Visualization

mean_percent_risk_in_gain_fair = mean([unlucky.percent_risk_in_gain_fair]);
mean_percent_risk_in_loss_fair = mean([unlucky.percent_risk_in_loss_fair]);
mean_percent_risk_in_mixed_fair = mean([unlucky.percent_risk_in_mixed_fair]);

mean_percent_risk_in_gain_biased = mean([unlucky.percent_risk_in_gain_biased]);
mean_percent_risk_in_loss_biased = mean([unlucky.percent_risk_in_loss_biased]);
mean_percent_risk_in_mixed_biased = mean([unlucky.percent_risk_in_mixed_biased]);

sem_percent_risk_in_gain_fair = sem([unlucky.percent_risk_in_gain_fair]);
sem_percent_risk_in_loss_fair = sem([unlucky.percent_risk_in_loss_fair]);
sem_percent_risk_in_mixed_fair = sem([unlucky.percent_risk_in_mixed_fair]);

sem_percent_risk_in_gain_biased = sem([unlucky.percent_risk_in_gain_biased]);
sem_percent_risk_in_loss_biased = sem([unlucky.percent_risk_in_loss_biased]);
sem_percent_risk_in_mixed_biased = sem([unlucky.percent_risk_in_mixed_biased]);

fair_percent_risk_means = [mean_percent_risk_in_gain_fair, mean_percent_risk_in_loss_fair, mean_percent_risk_in_mixed_fair]
biased_percent_risk_means = [mean_percent_risk_in_gain_biased,  mean_percent_risk_in_loss_biased, mean_percent_risk_in_mixed_biased]

fair_percent_risk_sems = [sem_percent_risk_in_gain_fair, sem_percent_risk_in_loss_fair, sem_percent_risk_in_mixed_fair];
biased_percent_risk_sems = [sem_percent_risk_in_gain_biased,  sem_percent_risk_in_loss_biased, sem_percent_risk_in_mixed_biased];

subplot(1,2,2)
labels = {'Gain','Loss','Mixed'};
b = bar([fair_percent_risk_means; biased_percent_risk_means]', 'grouped');
hold on
set(gca,'xtick',1:length(labels),'xticklabel',labels);
x = [b(1).XEndPoints; b(2).XEndPoints];
errorbar([b(1).XEndPoints; b(2).XEndPoints]', [fair_percent_risk_means; biased_percent_risk_means]', ...
        [fair_percent_risk_sems; biased_percent_risk_sems]', '.');
title(sprintf("Unlucky participants (n = %d)", N_UNLUCKY))
ylim([0 1])
ylabel('Percent risks taken');
legend({'Fair', 'Biased'})

p_lucky_gain = ranksum([lucky.percent_risk_in_gain_fair], [lucky.percent_risk_in_gain_biased]);
p_lucky_mixed = ranksum([lucky.percent_risk_in_mixed_fair], [lucky.percent_risk_in_mixed_biased]);
p_lucky_loss = ranksum([lucky.percent_risk_in_loss_fair], [lucky.percent_risk_in_loss_biased]);

p_unlucky_gain = ranksum([unlucky.percent_risk_in_gain_fair], [unlucky.percent_risk_in_gain_biased]);
p_unlucky_mixed = ranksum([unlucky.percent_risk_in_mixed_fair], [unlucky.percent_risk_in_mixed_biased]);
p_unlucky_loss = ranksum([unlucky.percent_risk_in_loss_fair], [unlucky.percent_risk_in_loss_biased]);

% CONTINUOUS

for i = 1 : N

    individual_data = data(i).data;

    gain_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) == 0);
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    data(i).percent_risk_in_gain = num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(:,4) == 0 & individual_data(:,5) < 0);
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    data(i).percent_risk_in_loss = num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) < 0);
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    data(i).percent_risk_in_mixed = num_mixed_gambles / length(mixed_trials);

end

data_table = struct2table(data); 
table_sorted = sortrows(data_table, 'percent_won_gambles'); 
data_sorted = table2struct(table_sorted);

luck = zeros(N, 1);
risk_in_gain = zeros(N, 1);
risk_in_loss = zeros(N, 1);
risk_in_mixed = zeros(N, 1);

for i = 1 : N

    individual_struct = data_sorted(i);

    luck(i) = individual_struct.percent_won_gambles;
    risk_in_gain(i) = individual_struct.percent_risk_in_gain;
    risk_in_loss(i) = individual_struct.percent_risk_in_loss;
    risk_in_mixed(i) = individual_struct.percent_risk_in_mixed;

end

figure
plot(luck,smooth(risk_in_gain))
axis square;
hold on; plot(luck,smooth(risk_in_loss))
hold on; plot(luck,smooth(risk_in_mixed))
xlabel("Percent won gambles")
ylabel("Percent risks taken")
ylim([0 1])
title(sprintf("All participants (n = %d)", N))
legend({'Gain', 'Loss', 'Mixed'})
Effect of luck on attractiveness of risky gamble
As gambles become more attractive, it follows that an individual will be more inclined toward risk-seeking behavior. Let's next investigate whether perceived luck pushes participants toward more risk-seeking or risk-averse behavior.
tgaincutoff = [1.5 1.75 1.9 2.1 2.4   2.6 3 3.2 4 4.5 5.1]; %riskygain/safe ratio
tmixedcutoff = [0.4 0.6 0.8 0.9 1.1   1.3 1.6 2 2.8 3.5 5.1]; %-riskygain/riskyloss ratio
tlosscutoff = tgaincutoff; %riskyloss/safe ratio

tgaincutoff = tgaincutoff(1:2:end);
tmixedcutoff = tmixedcutoff(1:2:end);
tlosscutoff = tgaincutoff;

L = length(tgaincutoff)-1;

% LUCKY

for s=1:N_LUCKY
    t=lucky(s).data(1:T/2,1:8);

    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_lucky=zeros(1,L); tmixedbin_lucky=zeros(1,L); tlossbin_lucky=zeros(1,L);

    for n=1:L

        lucky(s).tgainbin_fair(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        lucky(s).tmixedbin_fair(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        lucky(s).tlossbin_fair(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

    end

    t=lucky(s).data(T/2+1:T,1:8);

    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_lucky=zeros(1,L); tmixedbin_lucky=zeros(1,L); tlossbin_lucky=zeros(1,L);

    for n=1:length(tgaincutoff)-1

        lucky(s).tgainbin_biased(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        lucky(s).tmixedbin_biased(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        lucky(s).tlossbin_biased(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

    end

end

tgainbin_fair=vertcat(lucky.tgainbin_fair);
tmixedbin_fair=vertcat(lucky.tmixedbin_fair);
tlossbin_fair=vertcat(lucky.tlossbin_fair);

tgainbin_biased=vertcat(lucky.tgainbin_biased);
tmixedbin_biased=vertcat(lucky.tmixedbin_biased);
tlossbin_biased=vertcat(lucky.tlossbin_biased);

figure; 
subplot(2,3,1); errorbar(1:L,nanmean(tgainbin_fair),nanstd(tgainbin_fair)/sqrt(N_LUCKY)); axis square; 
hold on; errorbar(1:L,nanmean(tgainbin_biased),nanstd(tgainbin_biased)/sqrt(N_LUCKY)); ylabel("% gain gambles chosen"); legend({'Fair','Biased'},'location','southeast')
hold on
subplot(2,3,2); errorbar(1:L,nanmean(tlossbin_fair),nanstd(tlossbin_fair)/sqrt(N_LUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; errorbar(1:L,nanmean(tlossbin_biased),nanstd(tlossbin_biased)/sqrt(N_LUCKY)); ylabel("% loss gambles chosen"); 
hold on
subplot(2,3,3); errorbar(1:L,nanmean(tmixedbin_fair),nanstd(tmixedbin_fair)/sqrt(N_LUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; errorbar(1:L,nanmean(tmixedbin_biased),nanstd(tmixedbin_biased)/sqrt(N_LUCKY)); ylabel("% mixed gambles chosen");
hold on
for n=1:3, subplot(2,3,n); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Lucky participants (n = %d)',N_LUCKY)); end

% UNLUCKY

for s=1:N_UNLUCKY
    t=unlucky(s).data(1:T/2,1:8);

    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_lucky=zeros(1,L); tmixedbin_lucky=zeros(1,L); tlossbin_lucky=zeros(1,L);

    for n=1:length(tgaincutoff)-1

        unlucky(s).tgainbin_fair(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        unlucky(s).tmixedbin_fair(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        unlucky(s).tlossbin_fair(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

    end

    t=unlucky(s).data(T/2+1:T,1:8);

    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_lucky=zeros(1,L); tmixedbin_lucky=zeros(1,L); tlossbin_lucky=zeros(1,L);

    for n=1:length(tgaincutoff)-1

        unlucky(s).tgainbin_biased(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        unlucky(s).tmixedbin_biased(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        unlucky(s).tlossbin_biased(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

    end

end

tgainbin_fair=vertcat(unlucky.tgainbin_fair);
tmixedbin_fair=vertcat(unlucky.tmixedbin_fair);
tlossbin_fair=vertcat(unlucky.tlossbin_fair);

tgainbin_biased=vertcat(unlucky.tgainbin_biased);
tmixedbin_biased=vertcat(unlucky.tmixedbin_biased);
tlossbin_biased=vertcat(unlucky.tlossbin_biased);

subplot(2,3,4); errorbar(1:L,nanmean(tgainbin_fair),nanstd(tgainbin_fair)/sqrt(N_UNLUCKY)); axis square; 
hold on; errorbar(1:L,nanmean(tgainbin_biased),nanstd(tgainbin_biased)/sqrt(N_UNLUCKY)); ylabel("% gain gambles chosen"); legend({'Fair','Biased'},'location','southeast')
hold on
subplot(2,3,5); errorbar(1:L,nanmean(tlossbin_fair),nanstd(tlossbin_fair)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; errorbar(1:L,nanmean(tlossbin_biased),nanstd(tlossbin_biased)/sqrt(N_UNLUCKY)); ylabel("% loss gambles chosen"); 
hold on
subplot(2,3,6); errorbar(1:L,nanmean(tmixedbin_fair),nanstd(tmixedbin_fair)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; errorbar(1:L,nanmean(tmixedbin_biased),nanstd(tmixedbin_biased)/sqrt(N_UNLUCKY)); ylabel("% mixed gambles chosen");
hold on
for n=1:3, subplot(2,3,n+3); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Unlucky participants (n = %d)',N_UNLUCKY)); end
Model-based analysis
% LUCKY

for s=1:N_LUCKY
    t=lucky(s).data(T/2+1:T,1:11);

    result = fitmodel_prospect_theory(t); %fits 4-parameter prospect theory model
    lucky(s).result_pt = result;

    result = fitmodel_loss_aversion(t); %fits 2-parameter loss aversion model
    lucky(s).result_la = result;

    result = fitmodel_aa_model(t); %fits 6-parameter approach-avoidance model
    lucky(s).result_aa = result;

    result = fitmodel_prospect_theory_luck(t); %fits 5-parameter prospect theory model
    lucky(s).result_ptl = result;

    lucky(s).bic_pt = -2 * lucky(s).result_pt.modelLL + length(lucky(s).result_pt.b) * log(T);
    lucky(s).bic_la = -2 * lucky(s).result_la.modelLL + length(lucky(s).result_la.b) * log(T);
    lucky(s).bic_aa = -2 * lucky(s).result_aa.modelLL + length(lucky(s).result_aa.b) * log(T);
    lucky(s).bic_ptl = -2 * lucky(s).result_ptl.modelLL + length(lucky(s).result_ptl.b) * log(T);
    
    t(:,8) = lucky(s).result_pt.probchoice;
    t(:,9) = lucky(s).result_la.probchoice;
    t(:,10) = lucky(s).result_aa.probchoice;
    t(:,11) = lucky(s).result_ptl.probchoice;

    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_lucky=zeros(1,L); tmixedbin_lucky=zeros(1,L); tlossbin_lucky=zeros(1,L);

    for n=1:length(tgaincutoff)-1
        lucky(s).tgainbin(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        lucky(s).tmixedbin(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        lucky(s).tlossbin(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

        lucky(s).tgainbin_pt(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),8));
        lucky(s).tmixedbin_pt(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),8));
        lucky(s).tlossbin_pt(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),8));

        lucky(s).tgainbin_la(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),9));
        lucky(s).tmixedbin_la(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),9));
        lucky(s).tlossbin_la(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),9));

        lucky(s).tgainbin_aa(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),10));
        lucky(s).tmixedbin_aa(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),10));
        lucky(s).tlossbin_aa(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),10));

        lucky(s).tgainbin_ptl(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),11));
        lucky(s).tmixedbin_ptl(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),11));
        lucky(s).tlossbin_ptl(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),11));
    end

end

tgainbin_lucky=vertcat(lucky.tgainbin);
tmixedbin_lucky=vertcat(lucky.tmixedbin);
tlossbin_lucky=vertcat(lucky.tlossbin);

tgainbin_pt_lucky=vertcat(lucky.tgainbin_pt);
tmixedbin_pt_lucky=vertcat(lucky.tmixedbin_pt);
tlossbin_pt_lucky=vertcat(lucky.tlossbin_pt);

tgainbin_la=vertcat(lucky.tgainbin_la);
tmixedbin_la=vertcat(lucky.tmixedbin_la);
tlossbin_la=vertcat(lucky.tlossbin_la);

tgainbin_aa=vertcat(lucky.tgainbin_aa);
tmixedbin_aa=vertcat(lucky.tmixedbin_aa);
tlossbin_aa=vertcat(lucky.tlossbin_aa);

tgainbin_ptl=vertcat(lucky.tgainbin_ptl);
tmixedbin_ptl=vertcat(lucky.tmixedbin_ptl);
tlossbin_ptl=vertcat(lucky.tlossbin_ptl);

figure; 
subplot(2,3,1); errorbar(1:L,nanmean(tgainbin_lucky),nanstd(tgainbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction gain gambles chosen'); 
hold on; plot(1:L, nanmean(tgainbin_pt_lucky))
hold on; plot(1:L, nanmean(tgainbin_la))
hold on; plot(1:L, nanmean(tgainbin_aa))
legend({'Data','PT','LA','AA'},'location','southeast')

subplot(2,3,2); errorbar(1:L,nanmean(tlossbin_lucky),nanstd(tlossbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; plot(1:L, nanmean(tlossbin_pt_lucky))
hold on; plot(1:L, nanmean(tlossbin_la))
hold on; plot(1:L, nanmean(tlossbin_aa))

subplot(2,3,3); errorbar(1:L,nanmean(tmixedbin_lucky),nanstd(tmixedbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; plot(1:L, nanmean(tmixedbin_pt_lucky))
hold on; plot(1:L, nanmean(tmixedbin_la))
hold on; plot(1:L, nanmean(tmixedbin_aa))

for n=1:3, subplot(2,3,n); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Lucky participants (n = %d)',N_LUCKY)); end

% UNLUCKY

for s=1:N_UNLUCKY
    t=unlucky(s).data(T/2+1:T,1:11);

    result = fitmodel_prospect_theory(t); %fits 4-parameter prospect theory model
    unlucky(s).result_pt = result;

    result = fitmodel_loss_aversion(t); %fits 2-parameter loss aversion model
    unlucky(s).result_la = result;

    result = fitmodel_aa_model(t); %fits 6-parameter approach-avoidance model
    unlucky(s).result_aa = result;

    result = fitmodel_prospect_theory_luck(t); %fits 5-parameter prospect theory model
    unlucky(s).result_ptl = result;

    unlucky(s).bic_pt = -2 * unlucky(s).result_pt.modelLL + length(unlucky(s).result_pt.b) * log(T);
    unlucky(s).bic_la = -2 * unlucky(s).result_la.modelLL + length(unlucky(s).result_la.b) * log(T);
    unlucky(s).bic_aa = -2 * unlucky(s).result_aa.modelLL + length(unlucky(s).result_aa.b) * log(T);
    unlucky(s).bic_ptl = -2 * unlucky(s).result_ptl.modelLL + length(unlucky(s).result_ptl.b) * log(T);
    
    t(:,8) = unlucky(s).result_pt.probchoice;
    t(:,9) = unlucky(s).result_la.probchoice;
    t(:,10) = unlucky(s).result_aa.probchoice;
    t(:,11) = unlucky(s).result_ptl.probchoice;
    tgain = t(t(:,3)>0,:); tgain(:,6) = tgain(:,4)./tgain(:,3);
    tmixed = t(t(:,3)==0,:); tmixed(:,6) = -tmixed(:,4)./tmixed(:,5);
    tloss = t(t(:,3)<0,:); tloss(:,6) = tloss(:,5)./tloss(:,3);
    tgainbin_unlucky=zeros(1,L); tmixedbin_unlucky=zeros(1,L); tlossbin_unlucky=zeros(1,L);
    for n=1:length(tgaincutoff)-1
        unlucky(s).tgainbin(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),7));
        unlucky(s).tmixedbin(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),7));
        unlucky(s).tlossbin(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),7));

        unlucky(s).tgainbin_pt(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),8));
        unlucky(s).tmixedbin_pt(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),8));
        unlucky(s).tlossbin_pt(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),8));

        unlucky(s).tgainbin_la(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),9));
        unlucky(s).tmixedbin_la(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),9));
        unlucky(s).tlossbin_la(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),9));

        unlucky(s).tgainbin_aa(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),10));
        unlucky(s).tmixedbin_aa(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),10));
        unlucky(s).tlossbin_aa(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),10));

        unlucky(s).tgainbin_ptl(n)=mean(tgain(tgain(:,6)>tgaincutoff(n)&tgain(:,6)<tgaincutoff(n+1),11));
        unlucky(s).tmixedbin_ptl(n)=mean(tmixed(tmixed(:,6)>tmixedcutoff(n)&tmixed(:,6)<tmixedcutoff(n+1),11));
        unlucky(s).tlossbin_ptl(n)=mean(tloss(tloss(:,6)>tlosscutoff(n)&tloss(:,6)<tlosscutoff(n+1),11));
    end
end
tgainbin_unlucky=vertcat(unlucky.tgainbin);
tmixedbin_unlucky=vertcat(unlucky.tmixedbin);
tlossbin_unlucky=vertcat(unlucky.tlossbin);

tgainbin_pt_unlucky=vertcat(unlucky.tgainbin_pt);
tmixedbin_pt_unlucky=vertcat(unlucky.tmixedbin_pt);
tlossbin_pt_unlucky=vertcat(unlucky.tlossbin_pt);

tgainbin_la=vertcat(unlucky.tgainbin_la);
tmixedbin_la=vertcat(unlucky.tmixedbin_la);
tlossbin_la=vertcat(unlucky.tlossbin_la);

tgainbin_aa=vertcat(unlucky.tgainbin_aa);
tmixedbin_aa=vertcat(unlucky.tmixedbin_aa);
tlossbin_aa=vertcat(unlucky.tlossbin_aa);

subplot(2,3,4); errorbar(1:L,nanmean(tgainbin_unlucky),nanstd(tgainbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction gain gambles chosen');
hold on; plot(1:L, nanmean(tgainbin_pt_unlucky))
hold on; plot(1:L, nanmean(tgainbin_la))
hold on; plot(1:L, nanmean(tgainbin_aa))

subplot(2,3,5); errorbar(1:L,nanmean(tlossbin_unlucky),nanstd(tlossbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; plot(1:L, nanmean(tlossbin_pt_unlucky))
hold on; plot(1:L, nanmean(tlossbin_la))
hold on; plot(1:L, nanmean(tlossbin_aa))

subplot(2,3,6); errorbar(1:L,nanmean(tmixedbin_unlucky),nanstd(tmixedbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; plot(1:L, nanmean(tmixedbin_pt_unlucky))
hold on; plot(1:L, nanmean(tmixedbin_la))
hold on; plot(1:L, nanmean(tmixedbin_aa))

for n=1:3, subplot(2,3,n+3); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Unlucky participants (n = %d)',N_UNLUCKY)); end

figure; 
subplot(2,3,1); errorbar(1:L,nanmean(tgainbin_lucky),nanstd(tgainbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction gain gambles chosen'); 
hold on; plot(1:L, nanmean(tgainbin_pt_lucky))
hold on; plot(1:L, nanmean(vertcat(lucky.tgainbin_ptl)))
legend({'Data','PT','PTL'},'location','southeast')

subplot(2,3,2); errorbar(1:L,nanmean(tlossbin_lucky),nanstd(tlossbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; plot(1:L, nanmean(tlossbin_pt_lucky))
hold on; plot(1:L, nanmean(vertcat(lucky.tlossbin_ptl)))

subplot(2,3,3); errorbar(1:L,nanmean(tmixedbin_lucky),nanstd(tmixedbin_lucky)/sqrt(N_LUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; plot(1:L, nanmean(tmixedbin_pt_lucky))
hold on; plot(1:L, nanmean(vertcat(lucky.tmixedbin_ptl)))

for n=1:3, subplot(2,3,n); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Lucky participants (n = %d)',N_LUCKY)); end

subplot(2,3,4); errorbar(1:L,nanmean(tgainbin_unlucky),nanstd(tgainbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction gain gambles chosen');
hold on; plot(1:L, nanmean(tgainbin_pt_unlucky))
hold on; plot(1:L, nanmean(vertcat(unlucky.tgainbin_ptl)))

subplot(2,3,5); errorbar(1:L,nanmean(tlossbin_unlucky),nanstd(tlossbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction loss gambles chosen');
hold on; plot(1:L, nanmean(tlossbin_pt_unlucky))
hold on; plot(1:L, nanmean(vertcat(unlucky.tlossbin_ptl)))

subplot(2,3,6); errorbar(1:L,nanmean(tmixedbin_unlucky),nanstd(tmixedbin_unlucky)/sqrt(N_UNLUCKY)); axis square; ylabel('Fraction mixed gambles chosen');
hold on; plot(1:L, nanmean(tmixedbin_pt_unlucky))
hold on; plot(1:L, nanmean(vertcat(unlucky.tmixedbin_ptl)))

for n=1:3, subplot(2,3,n+3); axis([0 L+1 0 1]); xlabel('Gamble value quintile'); title(sprintf('Unlucky participants (n = %d)',N_UNLUCKY)); end

lucky_bic_pt = sum([lucky.bic_pt])
lucky_bic_la = sum([lucky.bic_la])
lucky_bic_aa = sum([lucky.bic_aa])
lucky_bic_ptl = sum([lucky.bic_ptl])


unlucky_bic_pt = sum([unlucky.bic_pt])
unlucky_bic_la = sum([unlucky.bic_la])
unlucky_bic_aa = sum([unlucky.bic_aa])
unlucky_bic_ptl = sum([unlucky.bic_ptl])


function s = sem(M)
    s = nanstd(M) / (sqrt(length(M)));
end

function b = get_bin(num,bins)
    for i = 2:length(bins)
        if num < bins(i)
            b = bins(i - 1);
            return
        end
    end
    b = -1;
end

function s = smooth(data)
    N_SMOOTH = 20;
    s = filter(ones(N_SMOOTH,1)./N_SMOOTH, 1, data);
end

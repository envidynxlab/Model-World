%% Distorion Index
%
% Model simulates washover deposition into an idealised built environment,
% for a range of built fractions.
%
% This model is an adaptation of the model detailed in:
%
% Lazarus ED, *Armstrong S (2015) Self-organized pattern formation in
% coastal barrier washover deposits, Geology 43(4), 363?366.
%
% This work is a contribution to a project supported by the Leverhulme
% Trust (RPG-2018-282) ? a collaboration with Evan Goldstein, Hannah
% Williams & Luke Taylor.
%
% Eli Lazarus (e.d.lazarus@soton.ac.uk)
% All rights reserved.
%
%%

clear all
close all
rand('seed', 0)

%% Setting up the domain:

N = 100; % edge length of square domain N^2
T = 20; % number of times the washover routine will iterate for a single run

DT = 1; % total time period (for diffusion step)
dt = 0.01; % time interval (for diffusion step)
K = 0.1; % diffusivity coefficient

rem = 0.1; % proportion of water discharge that gets eroded from bed (QS)
S = rem; % sediment lag (as portion of water flux)

P = 0.3; % default percent ht of berm that fails in a new failure

M = 10;
n = 2*M; % radius of influence for the Theim equation (relevant if multiple breaches)

Qo = 1; % initial water level (set-up against barrier)
level = Qo; % initial water level (set-up against barrier)
levelo = level; % store initial level
site = N/2; % default site for initial breach ? or see 'series' of sites, below

urbmax = 2; % elevation of non-erodible built environment


rrr = 1; % row counter, for variable storage

for spac_r = 4:15 % row spacing determines extent of built fraction
    
    spac_c = spac_r; % row spacing equals column spacing (square blocks of built environment)
    disp(['Spacer = ' num2str(spac_r)]) % display which step of the loop is active
    
    ccc = 1; % column counter, for variable storage
    for series = 1:25 % punches through the barrier in a bunch of different places...
        site = round(20 + (80 - 20).*rand); %...randomly selected within the middle 60 cells of barrier
        
        %% Puts buildings in the domain:
        
        Zurb = urbmax.*ones(N, N); % will become the built domain
        Z = zeros(N, N); % flat underlying domain
        
        for i = 0:3 % make streets 4 cells wide, slice up Zurb...
            Zurb((2 + i):spac_r:N, :) = 0;
            Zurb(:, (2 + i):spac_c:N) = 0;
        end
        
        Z = Zurb;
        
        %%
        
        Z(:, 1) = Qo; % make leading edge (first column, all rows) the barrier
        Zo = Z; % save original
        dZ = zeros(size(Z)); % domain to store temp changes in Z
        
        % Domain topography
        [col, row] = meshgrid(1:N, 1:N);
        
        berm = Z(:, 1); % isolates the barrier edge of the domain as its own alongshore profile
        
        
        if(1) % randomly determines breach depth for run as a proportion of barrier ht (bounded 0.1-0.7)
            
            P = 0.1 + (0.7 - 0.1).*rand;
            
        end
        
        berm(site) = 1 - P; % places the breach site in the berm
        bermo = berm; % store initial berm profile
        
        minwatero = 0.01; % min water threshold required to move sediment
        minwater = minwatero;
        
        totW = zeros(N, N); % stores cumulative total of cells "wet" by overwash flow
        WET = zeros(N, N); % wet cells
        WET2 = zeros(N, N); % temporary inventory of wet cells
        store_dips = zeros(N, T); % storage of any cells that "dip" to elevation < 0
        
        idx = 0; % counter
        %%
        for t = 1:T % for overwash "duration" of T steps
            
            tempberm = berm; % assigns berm height
            
            iter = 1; % counter
            for tt = dt:dt:DT % minor diffusion of the berm to knock off sharp corners ? equivalent to sand slumping down from edges; higher K increases effect
                
                del_z = dt*(K*(circshift(tempberm, 1) - 2*tempberm + circshift(tempberm, -1)));
                
                newz = tempberm + del_z;
                
                tempberm = newz;
                
                iter = iter + 1;
                
            end
            
            store_berm(:, t) = tempberm; % store what the berm now looks like
            store_level(t) = level; % store water level for run step
            
            Z(:, 1) = tempberm; % assign berm topo to barrier edge of domain
            
            
            %% Adjust water level according to berm changes
            % This whole section is book-keeping related to re-apportioning
            % overwash discharge through the breach(es)...
            
            sum_gap = sum(Qo - tempberm)/N; % how big is the breach?
            level = Qo - sum_gap; % adjust the water level across whole barrier accordingly
            
            %%
            qw = level - tempberm; % how much water moves through the breach
            qw(qw < 0) = 0; % no negatives
            
            [pks, locs] = findpeaks(qw); % find location of max discharge
            
            segs = zeros(size(tempberm)); % relevant breached segment(s) of the barrier
            prop = zeros(size(qw)); % for tracking water-discharge proportions
            
            if length(locs) == 1 % if there's only one breach site, redistribute proportion overwash discharge per cell
                
                port1 = locs/N;
                port2 = (N - locs)/N;
                
                segs(1:locs) = port1;
                segs((locs + 1):N) = port2;
                
                sum1 = sum(qw(1:locs));
                sum2 = sum(qw((locs + 1):N));
                
                prop(1:locs) = qw(1:locs)./sum1;
                prop((locs + 1):N) = qw((locs + 1):N)./sum2;
                
                rebalance = prop.*segs; % is the proportion of qw per cell
                
            else % if there are > 1 breach sites...
                
                % then water distribution goes from boundary to midpoint to
                % next breach OR midpoint to midpoint between two throats
                
                mps = locs(1:(end - 1)) + ceil(diff(locs)./2); % calc midpoints between breaches
                
                if ~isempty(mps)
                    for i = 1:(length(mps) + 1) % note: num segments == num mps + 1
                        
                        if i == 1
                            
                            port = mps(i)/N;
                            segs(1:mps) = port;
                            
                            sumi = sum(qw(1:mps(i)));
                            prop(1:mps(i)) = qw(1:mps(i))./sumi;
                            
                        elseif i == length(mps) + 1
                            
                            port = (N - mps(i - 1))/N;
                            segs((mps(i - 1) + 1):N) = port;
                            
                            sumi = sum(qw((mps(i - 1) + 1):N));
                            prop((mps(i - 1) + 1):N) = qw((mps(i - 1) + 1):N)./sumi;
                            
                        else
                            
                            port = (mps(i) - mps(i - 1))/N;
                            segs((mps(i - 1) + 1):mps(i)) = port;
                            
                            sumi = sum(qw((mps(i - 1) + 1):mps(i)));
                            prop((mps(i - 1) + 1):mps(i)) = qw((mps(i - 1) + 1):mps(i))./sumi;
                            
                        end
                        
                    end
                    
                    rebalance = prop.*segs; % is the proportion of qw per cell
                    
                end
                
            end
            
            rdwater = sum(qw).*rebalance; % redistributed water...
            store_rdw(:, t) = rdwater; % store variable
            
            
            %% Drawdown (cone of depression) sequence...
            % Irrelevant if there is only one breach site per run ?
            % otherwise, determines adjusted onshore forcing as a version
            % of a "drawdown" (Theim) problem...
            
            Q = level - tempberm;
            
            trans = 3; % in m/day
            A = 2*pi*trans; % Theim equation
            R = n; % radius of influence
            rw = 1:R;
            
            for i = 1:N
                
                shape = -Q(i)./A.*log(R./rw); % calculates for a radius out half the domain
                flip = fliplr(shape);
                
                if (i + n - 1) <= N
                    
                    sw(i, i:(i + n - 1)) = shape;
                    
                    if i >= n
                        
                        sw(i, (i - n + 1):i) = flip;
                        
                    else
                        
                        sw(i, 1:(i - 1)) = flip((n - i + 1):(n - 1));
                        
                    end
                    
                else
                    
                    sw(i, i:N) = shape(1:(N - i + 1));
                    
                    sw(i, (i - n + 1):i) = flip;
                    
                end
                
            end
            
            dd = sum(sw, 1); % sum all cones of depression (drawdown), get the profile
            
            
            %% Initiate overwash flow...
            
            Z(:, 1) = tempberm; % update barrier topography
            
            Ztemp = Z; % temporary variable, for adjustment
            Ztemp(:, 1) = bermo - rdwater; % adjust berm ht
            
            % Need to make sure the berm doesn't dig itself a pit...checks
            % the difference between the barrier & the first part of the
            % back-barrier domain
            dZ1Z2 = Z(:, 1) - Z(:, 2);
            dips = find(dZ1Z2 < 0);
            
            if ~isempty(dips)
                store_dips(dips, t) = 1;
            end
            
            %%
            Ztempo = Ztemp; % a couple of storage versions, for later checks
            Ztempoo = Ztemp;
            
            W = zeros(N, N); % difference between surge ht and Z ht; assume barrier toe == 0
            
            W(:, 1) = rdwater; % water ht
            W(W < 0) = 0; % no negatives in the water ht
            
            QA = W(:, 1); % initial discharge, stored as an array...
            
            move = zeros(N, N); % which cells are "moving" (water flow) at this step?
            move(W > 0) = 1; % just a black/white mask
            
            wet = zeros(N, N); % which cells are wet?
            wet(W > 0) = 1;
            
            % blanks for eventual surface update...
            add = zeros(size(Z));
            sub = zeros(size(Z));
            
            %%
            cycle = 0; % counter
            while sum(sum(move)) > 0 % if there are cells to move, will cycle until there aren't...
                
                nogo = find(Z >= level); % can't move into that cell
                move(nogo) = 0;
                
                cycle = cycle + 1; % update counter
                idx = idx + 1; % update counter
                
                Ztempo = Ztemp; % update "previous" terrain
                
                fronts = find(move > 0); % identify all the moving cells
                
                for f = 1:length(fronts) % for each moving cell...
                    
                    % where are you in the domain
                    r = row(fronts(f));
                    c = col(fronts(f));
                    
                    % water (if any) apportioned to each neigbour
                    qwdifset = zeros(8, 2); % water flux
                    qwdifset(:, 1) = 1:8;
                    
                    % sediment flux, as function of water, to each
                    % neighbour
                    qsdifset = zeros(8, 2); % sed flux
                    qsdifset(:, 1) = 1:8;
                    
                    % how to update the Z domain itself, based on sed flux qs
                    Zset = zeros(8, 2); % Z array
                    Zset(:, 1) = 1:8;
                    
                    
                    % Look around at neighbors (Queen's play [diagonals are game]), to determine if any
                    % water/sed goes there, and any update to Z domain...
                    
                    % upper L corner neighbour
                    if (r - 1) > 1 && (c - 1) > 1
                        
                        if wet(r - 1, c - 1) == 0
                            
                            qwdifset(1, 2) = (Ztemp(r, c) - Ztemp(r - 1, c - 1))/(sqrt(2));
                            qsdifset(1, 2) = (Z(r, c) - Z(r - 1, c - 1))/(sqrt(2));
                            Zset(1, 2) = Z(r - 1, c - 1);
                            
                        end
                        
                    end
                    
                    % neighbour immediately above
                    if (r - 1) > 1
                        
                        if wet(r - 1, c) == 0
                            
                            qwdifset(2, 2) = Ztemp(r, c) - Ztemp(r - 1, c);
                            qsdifset(2, 2) = Z(r, c) - Z(r - 1, c);
                            Zset(2, 2) = Z(r - 1, c);
                            
                        end
                        
                    end
                    
                    % upper R corner neighbour
                    if (r - 1) > 1 && (c + 1) < N
                        
                        if wet(r - 1, c + 1) == 0
                            
                            qwdifset(3, 2) = (Ztemp(r, c) - Ztemp(r - 1, c + 1))/(sqrt(2));
                            qsdifset(3, 2) = (Z(r, c) - Z(r - 1, c + 1))/(sqrt(2));
                            Zset(3, 2) = Z(r - 1, c + 1);
                            
                        end
                        
                    end
                    
                    % right neighbour
                    if (c + 1) < N
                        
                        if wet(r, c + 1) == 0
                            
                            qwdifset(4, 2) = Ztemp(r, c) - Ztemp(r, c + 1);
                            qsdifset(4, 2) = Z(r, c) - Z(r, c + 1);
                            Zset(4, 2) = Z(r, c + 1);
                            
                        end
                        
                    end
                    
                    
                    % lower right corner neighbour
                    if (r + 1) < N && (c + 1) < N
                        
                        if wet(r + 1, c + 1) == 0
                            
                            qwdifset(5, 2) = (Ztemp(r, c) - Ztemp(r + 1, c + 1))/(sqrt(2));
                            qsdifset(5, 2) = (Z(r, c) - Z(r + 1, c + 1))/(sqrt(2));
                            Zset(5, 2) = Z(r + 1, c + 1);
                            
                        end
                        
                    end
                    
                    % neighbour immediately below
                    if (r + 1) < N
                        
                        if wet(r + 1, c) == 0
                            
                            qwdifset(6, 2) = Ztemp(r, c) - Ztemp(r + 1, c);
                            qsdifset(6, 2) = Z(r, c) - Z(r + 1, c);
                            Zset(6, 2) = Z(r + 1, c);
                            
                        end
                        
                    end
                    
                    % lower left corner neighbour
                    if (r + 1) < N && (c - 1) > 1
                        
                        if wet(r + 1, c - 1) == 0
                            
                            qwdifset(7, 2) = (Ztemp(r, c) - Ztemp(r + 1, c - 1))/(sqrt(2));
                            qsdifset(7, 2) = (Z(r, c) - Z(r + 1, c - 1))/(sqrt(2));
                            Zset(7, 2) = Z(r + 1, c - 1);
                            
                        end
                        
                    end
                    
                    % left neighbour
                    if (c - 1) > 1
                        
                        if wet(r, c - 1) == 0
                            
                            qwdifset(8, 2) = Ztemp(r, c) - Ztemp(r, c - 1);
                            qsdifset(8, 2) = Z(r, c) - Z(r, c - 1);
                            Zset(8, 2) = Z(r, c - 1);
                            
                        end
                        
                    end
                    
                    
                    %%
                    qwdifseto = qwdifset; % preserve original (water flux)
                    qsdifseto = qsdifset; % preserve original (sed flux)
                    
                    qwpossum = 0; % intitalize positive water flux
                    for i = 1:8 % total positive difference (for proportional divvy)
                        
                        if qwdifset(i, 2) > 0 % keep positives
                            
                            qwpossum = qwpossum + qwdifset(i, 2);
                            
                        elseif qwdifset(i, 2) < 0 % clear negatives
                            
                            qwdifset(i, 2) = 0;
                            
                        end
                    end
                    
                    qspossum = 0; % initialize positive sediment flux
                    for i = 1:8 % total pos difference (for proportional divvy)
                        
                        if qsdifset(i, 2) > 0 % keep positives
                            
                            qspossum = qspossum + qsdifset(i, 2);
                            
                        elseif qsdifset(i, 2) < 0 % clear negatives
                            
                            qsdifset(i, 2) = 0;
                            
                        end
                    end
                    
                    qwsorted = sortrows(qwdifset, -2); % sort second column (descending)
                    qssorted = sortrows(qsdifset, -2); % sort second column (descending)
                    Zsorted = sortrows(Zset, -2); % Zset comes from Z (not Ztemp)
                    
                    % Need to move water ONLY down to the ht of min neighbor
                    qwvals = qwsorted(:, 2);
                    qwyes = find(qwvals > 0);
                    
                    Zvals = Zsorted(:, 2); %...likewise, comes from Z (not Ztemp)
                    
                    if ~isempty(qwyes) % 'not empty' means there ARE positive values for distribution...
                        
                        qwneighs = qwvals(qwyes);
                        qwn_min = min(qwneighs);
                        
                        if Z(r, c) >= min(Zvals) % empty all the water
                            
                            qwpossumo = qwpossum;
                            
                            perc = 1;
                            
                            QW = qwsorted;
                            QW(1:length(qwneighs), 2) = qwneighs;
                            
                        else
                            
                            qwpossumo = qwpossum;
                            
                            qwpossum = qwpossumo - qwn_min;
                            
                            perc = qwpossum/qwpossumo;
                            
                            % Now, drop the minimum from consideration...
                            
                            qwneighs(qwneighs == qwn_min) = 0;
                            
                            QW = qwsorted;
                            QW(1:length(qwneighs), 2) = qwneighs;
                            
                        end
                        
                    end
                    
                    QW = qwsorted;
                    qwdistrib = zeros(size(QW));
                    qwdistrib(:, 1) = QW(:, 1);
                    qsdistrib = qssorted;
                    
                    if qwpossum > minwater % IF there's (enough) water to move around...
                        
                        tempW = W;
                        
                        qwdistrib(:, 2) = W(r, c).*(QW(:, 2)./qwpossum);
                        
                        qsdistrib(:, 2) = S.*qwdistrib(:, 2); % sed lag from qw
                        
                        
                        % Don't dig a pit!
                        if c == 1
                            if ~isempty(dips)
                                
                                yes = find(r == dips);
                                
                                if ~isempty(yes)
                                    qsdistrib(:, 2) = 0;
                                end
                            end
                        end
                        
                        % Step through these inventories of neighbours...if there's
                        % enough water (and sed) to move, then move it to
                        % the appropriate cell. Move down the list
                        % according to rank...
                        for ii = 1:8
                            
                            if qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 1
                                
                                W(r - 1, c - 1) = W(r - 1, c - 1) + qwdistrib(ii, 2);
                                
                                add(r - 1, c - 1) = add(r - 1, c - 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r - 1, c - 1) = 1; % update the 'moving' domain
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 2
                                
                                W(r - 1, c) = W(r - 1, c) + qwdistrib(ii, 2);
                                
                                
                                add(r - 1, c) = add(r - 1, c) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r - 1, c) = 1;
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 3
                                
                                W(r - 1, c + 1) = W(r - 1, c + 1) + qwdistrib(ii, 2);
                                
                                
                                add(r - 1, c + 1) = add(r - 1, c + 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r - 1, c + 1) = 1;
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 4
                                
                                W(r, c + 1) = W(r, c + 1) + qwdistrib(ii, 2);
                                
                                
                                add(r, c + 1) = add(r, c + 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r, c + 1) = 1;
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 5
                                
                                W(r + 1, c + 1) = W(r + 1, c + 1) + qwdistrib(ii, 2);
                                
                                add(r + 1, c + 1) = add(r + 1, c + 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r + 1, c + 1) = 1;
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 6
                                
                                W(r + 1, c) = W(r + 1, c) + qwdistrib(ii, 2);
                                
                                add(r + 1, c) = add(r + 1, c) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                move(r + 1, c) = 1;
                                
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 7
                                
                                W(r + 1, c - 1) = W(r + 1, c - 1) + qwdistrib(ii, 2);
                                
                                
                                add(r + 1, c - 1) = add(r + 1, c - 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r + 1, c - 1) = 1;
                                
                            elseif qwdistrib(ii, 2) > 0 && qwdistrib(ii, 1) == 8
                                
                                W(r, c - 1) = W(r, c - 1) + qwdistrib(ii, 2);
                                
                                add(r, c - 1) = add(r, c - 1) + qsdistrib(ii, 2);
                                sub(r, c) = add(r, c) - qsdistrib(ii, 2);
                                
                                
                                move(r, c - 1) = 1;
                                
                            end
                            
                            
                            
                        end
                        
                        W(r, c) = 0; % empty the water cell
                        move(r, c) = 0; % reset cell that moved...
                        
                        wet(r, c) = 1; % % capture 'wet' cell
                        
                        
                        dZ = Z - Zo; % update overall domain
                        
                        
                    end % end the loop of any mover
                end % end loop of ALL movers in domain
                
                totW = totW + W; % update cumulative total 'wet' across domain
                
                Ztemp = W + Z; % water as an elevated surface on domain...
                Ztemp(:, 1) = Ztempoo(:, 1); % capture the intermediate domain
                
                WET = WET + wet; % update 'wet' domain
                
                
                if max(max(Ztemp - Ztempo)) == 0 % kick out of loop if there's no change
                    break
                end
                
            end % closes 'while' loop
            
            store_cycle(1, t) = cycle;
            
            
            WET2 = WET2 + wet; % updates OVERALL 'wet' domain
            
            net = add + sub; % map net change
            
            net(move == 1) = 0; % where 'move' cells landed, there's no change
            
            if ~isempty(dips)
                net(dips) = 0; % fix any dips
            end
            
            sn = net;
            sn(sn ~= 0) = t; % captures the cells that changed
            
            Z = Z + net; % update domain
            Z(Z < 0) = 0; % again, no holes...
            Z(:, (end - 5):end) = 0; % clear the end of the domain (though nothing should get close)
            
            
            store_sn(:, :, t) = sn; % stacks the net changers
            
            store_Z(:, :, t) = Z; % stacks the domain
            
            
            
            
            %% Make berm failure a 50/50 chance each t step
            if(0) % only relevant for multiple simultaneous breaches...
                nextfaild = 0;
                coin = round(rand);
                store_coin(t) = coin; % and track the sequence...
                
                if(coin)
                    
                    [ddpks, ddlocs] = findpeaks(dd);
                    
                    if ~isempty(ddpks)
                        
                        nfi = find(ddpks == max(ddpks));
                        
                        rp = randperm(length(nfi));
                        
                        nextfail = ddlocs(nfi(rp(1)));
                        
                        nextfaild = P*(tempberm(nextfail) - Z(nextfail, 2));
                        
                        if nextfaild < 0
                            nextfaild = 0;
                        end
                        
                        store_nfi(t) = nextfail;
                        
                    else
                        
                        [ddpks, ddlocs] = findpeaks(tempberm);
                        
                        nfi = find(ddpks == max(ddpks));
                        
                        nextfail = ddlocs(nfi);
                        
                        nextfaild = P*(tempberm(nextfail) - Z(nextfail, 2));
                        
                        
                        nextfail = 1;
                        nextfaild = 0;
                        
                        store_nfi(t) = nextfail;
                        
                        
                    end
                end
            end
            
            
            berm = Z(:, 1); % update barrier
            berm(berm < 0) = 0; % no negatives
            
            store_berm(:, t) = berm; % store alongshore barrier profile
            store_QA(:, t) = QA;
            
            store_move(:, :, t) = move;
            
            
        end % end of T step
        
        
        %% Intermediate ANALYTICS (tracked while model runs)
        
        Zdepo = Z;
        Zdepo = Zdepo - Zo; % map of all deposited material
        Zdepo(Zdepo < 0) = 0; % no negatives
        
        ZdepoBW = Zdepo;
        ZdepoBW(ZdepoBW > 0) = 5; % footprint of deposition, set high (> Zurb) to make it appear brighter than "buildings"
        ZDEPO = ZdepoBW + Zo; % put the change footprint atop intial domain
        
        ZdepoBW(ZdepoBW > 0) = 1; % depo footprint black/white...
        
        % Calculate intrusion length L:
        for i = 1:N
            
            if ~isempty(find(ZdepoBW(:, i) == 1))
                
                L = i;
                
            end
        end
        
        % Calculate (& store) other morphometric stats:
        
        STATS = regionprops(ZdepoBW, 'area', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter');
        area = STATS.Area;
        major = round(STATS.MajorAxisLength);
        minor = round(STATS.MinorAxisLength);
        perim = round(STATS.Perimeter);
        vol = round(sum(sum(Zdepo)));
        
        store_area(rrr, ccc) = area;
        store_perim(rrr, ccc) = perim;
        store_intrusion(rrr, ccc) = L;
        store_P(rrr, ccc) = P;
        
        Sden = length(find(Zo == 0))./(N^2);
        store_Sd(rrr, ccc) = Sden;
        
        SL = length(find(Zo == 0))/4;
        store_SL(rrr, ccc) = SL;
        
        Bden = length(find(Zo == 2))./(N^2);
        store_Bd(rrr, ccc) = Bden;
        
        
        ccc = ccc + 1; % iterate column index
                
    end % of the breach location 'series'
    
    rrr = rrr + 1; % iterate row index
    
end % of the urban fabric 'spacing' loop


%% DISTORTION INDEX

% Calculate distortion index (DI):
Pideal = (pi + 2).*sqrt(2.*store_area./pi);
DI = store_perim./Pideal;


%% FIGURES

figure
hold on
for i = 1:size(DI, 1)
    scatter(log(store_area(i, :)), DI(i, :), 50, (store_Bd(i, :)), 'filled')
end
hold off

hold on
for i = 1:5 % range covered by empirical data
    scatter(log(store_area(i, :)), DI(i, :), 50, 'k') 
end
hold off

box on
axis square
colorbar
ylabel('DI')
xlabel('Area')
title('URB DI v A by urban density (color)')


%%
figure
hold on
for i = 1:size(DI, 1)
    
    scatter(log(store_area(i, :)), log(store_perim(i, :)), 50, (store_Bd(i, :)), 'filled')
    
end

hold on
for i = 1:5 % range covered by empirical data
    scatter(log(store_area(i, :)), log(store_perim(i, :)), 50, 'k')
end
hold off
hold off
box on
axis square
colorbar
ylabel('log perimeter')
xlabel('log area')
title('model P vs A by BD (color)')



%%
pips = find(DI < 2.5); % pick out empirical range

figure
hold on
for i = 1:size(DI, 1)
    scatter(log(store_area(i, :)), log(store_perim(i, :)), 50, DI(i, :), 'filled')
        
end

scatter(log(store_area(pips)), log(store_perim(pips)), 50, 'k') % range covered by empirical data

hold off
box on
axis square
colorbar
ylabel('log perimeter')
xlabel('log area')
title('model P vs A by DI (color)')


%%
figure
hold on
plot(log(store_area(6:end, :)), log(store_intrusion(6:end, :)), '.k', 'MarkerSize', 15)
plot(log(store_area(2:5, :)), log(store_intrusion(2:5, :)), '.g', 'MarkerSize', 15) %BDs in real range
plot(log(store_area(1, :)), log(store_intrusion(1, :)), '.r', 'MarkerSize', 15) %BD = 0
hold off
box on
axis square
% colorbar
ylabel('log length')
xlabel('log area')


%%
figure
plot(store_Bd, DI, '.k', 'MarkerSize', 15)
box on
axis square
ylabel('DI')
xlabel('building density')


%%

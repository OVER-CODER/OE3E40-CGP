clc;
clear all;
%seq = load('human22_cleane_onlyACGT.mat');
%seq = seq.hnew;
seq = textread('AB016625.dat','%s');
seq = [seq{:}];
start_stop = cpgisland(seq);
starts_cpg = start_stop.Starts;
stops_cpg = start_stop.Stops;

for x = 1:3

    r=randperm(length(starts_cpg)); % to get random index numbers 
    train_starts_cpg = starts_cpg(r(1:0.9*length(starts_cpg))); % 90 percent cpg in training
    train_stops_cpg = stops_cpg(r(1:0.9*length(starts_cpg))); % 90 percent cpg in training

    test_starts_cpg = starts_cpg(r(0.9*length(starts_cpg)+1:length(r))); % remaining 10 percent cpg in test
    test_stops_cpg = stops_cpg(r(0.9*length(starts_cpg)+1:length(r))); % remaining 10 percent cpg in test



    init_pb = [];
    A_counts = 0;
    C_counts = 0;
    G_counts = 0;
    T_counts = 0;

    for i = 1:length(starts_cpg)

       if seq(starts_cpg(i))=='A'
           A_counts = A_counts+1;
       elseif seq(starts_cpg(i)) == 'C'
           C_counts = C_counts+1;
       elseif seq(starts_cpg(i)) == 'G'
           G_counts = G_counts+1;
       elseif seq(starts_cpg(i)) == 'T'
           T_counts = T_counts+1;

       end
    end


    init_pb(:,1) = A_counts./length(starts_cpg);
    init_pb(:,2) = C_counts./length(starts_cpg);
    init_pb(:,3) = G_counts./length(starts_cpg);
    init_pb(:,4) = T_counts./length(starts_cpg);



    cpg1 = seq(train_starts_cpg(1):train_stops_cpg(1)); % 1st random training cpg island start,stop values used for getting a cpg island

    dc1 = dimercount(cpg1);


    


    fields = fieldnames(dc1);  
    % transition_cpg1 will contain matrix
    for fieldname = fields'
        field = fieldname{1};
        value = dc1.(fieldname{1});

        if field(1) == 'A'
            total_possiblities_count = 0;
            total_possiblities_count = total_possiblities_count+dc1.AA;
            total_possiblities_count=total_possiblities_count+dc1.AC;
            total_possiblities_count=total_possiblities_count+dc1.AG;
            total_possiblities_count=total_possiblities_count+dc1.AT;

            if field(2) == 'A'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{1,1} = pb;
            elseif field(2) == 'C'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{1,2} = pb;
            elseif field(2) == 'G'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{1,3} = pb;
            elseif field(2) == 'T'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{1,4} = pb;
            end

        elseif field(1) == 'C'
            total_possiblities_count = 0;
            total_possiblities_count = total_possiblities_count+dc1.CA;
            total_possiblities_count=total_possiblities_count+dc1.CC;
            total_possiblities_count=total_possiblities_count+dc1.CG;
            total_possiblities_count=total_possiblities_count+dc1.CT;
            if field(2) == 'A'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{2,1} = pb;
            elseif field(2) == 'C'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{2,2} = pb;
            elseif field(2) == 'G'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{2,3} = pb;
            elseif field(2) == 'T'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{2,4} = pb;
            end

        elseif field(1) == 'G'
            total_possiblities_count = 0;
            total_possiblities_count = total_possiblities_count+dc1.GA;
            total_possiblities_count=total_possiblities_count+dc1.GC;
            total_possiblities_count=total_possiblities_count+dc1.GG;
            total_possiblities_count=total_possiblities_count+dc1.GT;
            if field(2) == 'A'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{3,1} = pb;
            elseif field(2) == 'C'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{3,2} = pb;
            elseif field(2) == 'G'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{3,3} = pb;
            elseif field(2) == 'T'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{3,4} = pb;
            end
        elseif field(1) == 'T'
            total_possiblities_count = 0;
            total_possiblities_count = total_possiblities_count+dc1.TA;
            total_possiblities_count=total_possiblities_count+dc1.TC;
            total_possiblities_count=total_possiblities_count+dc1.TG;
            total_possiblities_count=total_possiblities_count+dc1.TT;
            if field(2) == 'A'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{4,1} = pb;
            elseif field(2) == 'C'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{4,2} = pb;
            elseif field(2) == 'G'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{4,3} = pb;
            elseif field(2) == 'T'
                curr_dimer_count = value;
                pb = curr_dimer_count /  total_possiblities_count;
                transition_cpg1{4,4} = pb;
            end

        end
    end

    fprintf('Initial Probability for %d\n', x);
    disp(init_pb);
    fprintf('Transition Probability for %d\n', x);
    disp(transition_cpg1);
end





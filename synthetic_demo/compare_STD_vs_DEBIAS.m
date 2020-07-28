% this assumes that DEMO_20_06_synthetic_field.m has been run
addpath('../../3_SharedCode/')
assert(exist('DEMO_20_06_synthetic_field_COMPLETE','var') && ...
    DEMO_20_06_synthetic_field_COMPLETE, ...
    "Run DEMO_20_06_synthetic_field first")


figure(50); clf;
h =gcf;

Data_SV = svd(DataAssembled);

cmap = [winter(8);flipud(autumn(8)) ];
cmap_diff = [summer;flipud(autumn)];

METH1_NAME = "ROM DMD DEBIAS";
METH2_NAME = "ROM DMD EXACT";
Ns = [2:20, 20:5:80];

COMPARISON_NAME = genvarname(METH1_NAME) + "_VS_" + genvarname(METH2_NAME);
filename = "comp_" + COMPARISON_NAME + datestr(now,'yymmddHHMMSS');


ERROR_METH1 = nan(1,numel(Ns));
ERROR_METH2 = nan(1,numel(Ns));
ERROR_METH1_METH2 = nan(1,numel(Ns));


v = VideoWriter(filename,'MPEG-4');
v.FrameRate= 2;
v.Quality=95;
open(v);

REFERENCE = DataWoNoise;
REFERENCE_NAME = "Data w/o Noise";
NOISELEVEL = norm(NoiseComponent,'fro');

for k = 1:numel(Ns)

    %% METHOD1 - computation of modes and reconstruction
    out_METH1 = dmd(DataAssembled, dt, Ns(k),'debias',Ns(k) );
    [ROM_METH1,out_METH1] = reduce_order(out_METH1.Phi, out_METH1.omega, out_METH1.b,transpose(t), (1:Ns(k))  );

    %% METHOD2 - computation of modes and reconstruction
    out_METH2 = dmd(DataAssembled, dt, Ns(k) );
    [ROM_METH2,out_METH2] = reduce_order(out_METH2.Phi, out_METH2.omega, out_METH2.b,transpose(t), (1:Ns(k))   );

    % if reconstruction does not take into account conjugate modes, skip
    % this is easier than checking if two consecutive eigenvalues are
    % conjuage
    if  norm(imag(ROM_METH1)) > 1e-8 || ...
            norm(imag(ROM_METH2)) > 1e-8
        continue
    end
    
    % discard imaginary component - if the above condition was true, it's small
    % anyway
    ROM_METH1 = real(ROM_METH1);
    ROM_METH2 = real(ROM_METH2);

    %% Compute differences between models and models and data
    METH1_to_REF = REFERENCE - ROM_METH1;
    METH2_to_REF = REFERENCE - ROM_METH2;
    METH1_to_METH2 = ROM_METH1 - ROM_METH2;
    
    if k == 1
    clim_rec = max(abs(REFERENCE),[],'all') * [-1,1];
    clim_err = max( ...
        max(abs(METH1_to_REF),[],'all'),...
        max(abs(METH2_to_REF),[],'all')...        
         )* [-1,1];
    end

    %% Normalize summary errors
    Normalization = norm(REFERENCE,'fro');
    ERROR_METH1(k) = norm(METH1_to_REF,'fro')/Normalization;
    ERROR_METH2(k) = norm(METH2_to_REF,'fro')/Normalization;
    ERROR_METH1_METH2(k) = norm(METH1_to_METH2,'fro')/Normalization;


    %% VISUALIZATION

    h_true = subplot(3,3,1);
    pcolor(t, SpaceDomain, REFERENCE ); shading flat;
    colormap(cmap);
    colorbar;
    title("Reference: " + REFERENCE_NAME);
    caxis(h_true, clim_rec);

    h_meth1 = subplot(3,3,2);
    pcolor(t, SpaceDomain, ROM_METH1 ); shading flat;    
    colormap(cmap);
    colorbar;
    title({"METHOD1: " + METH1_NAME + " Model rank = " + out_METH1.model_rank  ;"Reconstruction rank = " + out_METH1.data_rank});
    caxis(h_meth1, clim_rec);

    h_error1 = subplot(3,3,3);
    h_DIFF = pcolor(t, SpaceDomain, METH1_to_REF ); shading flat;
    title({"METHOD1: " + METH1_NAME,"Difference to REF"})
    colorbar;
    colormap(gca, cmap_diff);
    caxis(h_error1, clim_err);

    h_meth2 = subplot(3,3,5);
    pcolor(t, SpaceDomain, ROM_METH2 ); shading flat;
    colorbar;
    title({"METHOD2: " + METH2_NAME + " Model rank = " + out_METH2.model_rank  ;"Reconstruction rank = " + out_METH2.data_rank});
    colormap(gca,cmap);
    caxis(h_meth2, clim_rec);


    h_error2 = subplot(3,3,6);
    pcolor(t, SpaceDomain, METH2_to_REF ); shading flat;
    title({"METHOD2: " + METH2_NAME,"Difference to REF"})
    colormap(gca, cmap_diff);
    colorbar;
    caxis(h_error2, clim_err);

    

    subplot(3,3,8);
    pcolor(t, SpaceDomain, METH1_to_METH2 ); shading flat;
    title("Difference between methods")
    colormap(gca, cmap_diff );
    colorbar;
    caxis([-1,1]*max(abs(METH1_to_METH2),[], 'all'));

    h_comp = subplot(3,3,9);

    plot(Ns(1:k), ERROR_METH1(1:k), 'r-o','DisplayName',"METH1 to REF"); hold on;
    plot(Ns(1:k), ERROR_METH2(1:k), 'b-x','DisplayName',"METH2 to REF" );
    plot(Ns(1:k), ERROR_METH1_METH2(1:k),'k-+', 'DisplayName','Between METH');
    hold off;
    xlim([0, max(Ns)]);
    ylim([0, max(1, NoiseMagnitudeGain*1.5)]);
    grid minor;
    set(h_comp,'yscale','log');

    legend('fontsize',6,'location','best');
    hy = yline(NOISELEVEL/Normalization,'r:', 'Noise Level');
    hy.Annotation.LegendInformation.IconDisplayStyle = 'off'; % skip in legend;

    title("Error as REF-normalized Frob. norm between field matrices");


    subplot(3,3,4);cla;
    scatter( real(E), imag(E), 100+10*randn(size(E)), 'ko', 'DisplayName','JAC' ); hold on;
    axis([-0.5,0.5, -7,7]);
    scatter( real(out_METH1.omega), imag(out_METH1.omega), 100, 'rx', 'DisplayName',"METH1" );
    scatter( real(out_METH2.omega), imag(out_METH2.omega), 100, 'b+', 'DisplayName',"METH2" ); hold off;
    legend('fontsize',6,'location','southwest');
    grid minor;

    title("Three spectra - zoomed in")
    xlabel("Re(\omega)")
    ylabel("Im(\omega)")

    subplot(3,3,7); cla;
    copyobj(get(subplot(3,3,4),'Children'),...
        subplot(3,3,7))     % copy spectrum plot verbatim
    grid minor;
    axis auto; % zoom out;
    title("Three spectra - zoomed out")
    xlabel("Re(\omega)")
    ylabel("Im(\omega)")

    sgtitle({"Reduced order models; Dim = " + Ns(k) + " TimeScale = " + TimeScale; ...
        IncludedComponents; ExcludedComponents;
        },'horizontalAlignment', 'left' );


    frame = getframe(h);
    writeVideo(v,frame);

end

close(v);

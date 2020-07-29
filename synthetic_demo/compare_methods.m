% this assumes that synthetic_field.m has been run
assert(exist('synthetic_field_COMPLETE','var') && ...
    synthetic_field_COMPLETE, ...
    "Run synthetic_field.m first")


figure(50); clf;
h =gcf;

Data_SV = svd(DataAssembled);

cmap = [winter(8);flipud(autumn(8)) ];
cmap_diff = [summer;flipud(autumn)];

METH1_NAME = "TDMD";
METH2_NAME = "TDMD Trunc Input";
Ns = [2:2:20, 20:5:50];
Ns = repmat(30, 1, 20);
idx = floor( linspace(10, size(DataAssembled,2), numel(Ns)) );
Ns = min(idx-1, Ns);

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
    
    
    METH1_fitstep = [1,idx(k)];
    METH2_fitstep = [1,idx(k)];
    
    
    %% METHOD1 - computation of modes and reconstruction
    out_METH1 = dmd(DataAssembled, dt, 95, 'rom_type', 'tlsq','step',METH1_fitstep );
    
    [ROM_METH1,out_METH1] = reduce_order(out_METH1.Phi, out_METH1.omega, out_METH1.b,transpose(t), (1:Ns(k))  );
    
    %% METHOD2 - computation of modes and reconstruction
    out_METH2 = dmd(DataAssembled(:,1:idx(k)), dt, Ns(k), 'rom_type', 'tlsq','step',METH2_fitstep );
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
    
    
    %% Normalize summary errors
    Normalization = norm(REFERENCE,'fro');
    ERROR_METH1(k) = norm(METH1_to_REF,'fro')/Normalization;
    ERROR_METH2(k) = norm(METH2_to_REF,'fro')/Normalization;
    ERROR_METH1_METH2(k) = norm(METH1_to_METH2,'fro')/Normalization;
    
    
    %% VISUALIZATION
    balance_color = @(hh,ff)caxis(hh, max(abs(ff),[],'all')*[-1,1]);
    
    
    h_true = subplot(3,3,1);
    pcolor(t, SpaceDomain, REFERENCE ); shading flat;
    colormap(cmap);
    colorbar;
    title("Reference: " + REFERENCE_NAME,'interpreter','none');
    balance_color(h_true, REFERENCE);
    
    h_meth1 = subplot(3,3,2);
    pcolor(t, SpaceDomain, ROM_METH1 ); shading flat;
    colormap(cmap);
    colorbar;
    title({"METHOD1: " + METH1_NAME + " Model rank = " + out_METH1.model_rank  ;"Reconstruction rank = " + out_METH1.data_rank},...
        'interpreter','none');
    balance_color(h_meth1, REFERENCE);
    for k = 1:numel(METH1_fitstep)
        hx1(k) = xline( (METH1_fitstep(k)-1)*dt, 'k:','Fit here');
        hx1(k).LineWidth = 3;
    end
    
    
    h_error1 = subplot(3,3,3);
    h_DIFF = pcolor(t, SpaceDomain, METH1_to_REF ); shading flat;
    title({"METHOD1: " + METH1_NAME,"Difference to REF"},'interpreter','none')
    colorbar;
    colormap(gca, cmap_diff);
    balance_color(h_error1, METH1_to_REF);
    
    h_meth2 = subplot(3,3,5);
    pcolor(t, SpaceDomain, ROM_METH2 ); shading flat;
    colorbar;
    title({"METHOD2: " + METH2_NAME + " Model rank = " + out_METH2.model_rank  ;"Reconstruction rank = " + out_METH2.data_rank},...
        'interpreter','none');
    colormap(gca,cmap);
    balance_color(h_meth2, REFERENCE);
    for k = 1:numel(METH2_fitstep)
        hx2(k) = xline( (METH2_fitstep(k)-1)*dt, 'k:','Fit here');
        hx2(k).LineWidth = 3;
    end
    
    
    h_error2 = subplot(3,3,6);
    pcolor(t, SpaceDomain, METH2_to_REF ); shading flat;
    title({"METHOD2: " + METH2_NAME,"Difference to REF"},'interpreter','none')
    colormap(gca, cmap_diff);
    colorbar;
    balance_color(h_error2, METH2_to_REF);
    
    
    
    h_error3 = subplot(3,3,8);
    pcolor(t, SpaceDomain, METH1_to_METH2 ); shading flat;
    title("Difference between methods")
    colormap(gca, cmap_diff );
    colorbar;
    balance_color(h_error3, METH1_to_METH2);
    
    h_comp = subplot(3,3,9); cla;
    
    h_nerr1 = plot(idx, ERROR_METH1, 'r-o','DisplayName',"METH1 to REF"); hold on;
    h_nerr2 = plot(idx, ERROR_METH2, 'b-x','DisplayName',"METH2 to REF" );
    h_nerr3 = plot(idx, ERROR_METH1_METH2,'k-+', 'DisplayName','Between METH');
    hold off;
    xlim([0, max(idx)]);
    grid minor;
    set(h_comp,'yscale','log');
    
    legend('fontsize',6,'location','best');
    hy = yline(NOISELEVEL/Normalization,'r:', 'Noise Level');
    hy.Annotation.LegendInformation.IconDisplayStyle = 'off'; % skip in legend;
    
    title("Normalized Frob. norm error");
    xlabel("Truncation index");
    
    
    subplot(3,3,4);cla;
    scatter( real(E), imag(E), 100+10*randn(size(E)), 'ko', 'DisplayName','JAC' ); hold on;
    ylim([-7,7]);
    xlims = xlim;
    xlims(1) = -log(2)/( max(TimeDomain)/10 );
    xlim(xlims);
    
    sizing = @(x)50*normalize(abs(x));
    scatter( real(out_METH1.omega), imag(out_METH1.omega), max( 100+sizing(out_METH1.b), 10), 'rx', 'DisplayName',"METH1" );
    scatter( real(out_METH2.omega), imag(out_METH2.omega), max( 100+sizing(out_METH2.b), 10), 'b+', 'DisplayName',"METH2" ); hold off;
    legend('fontsize',6,'location','best');
    grid minor;
    
    title("Three spectra - zoomed in")
    xlabel("Re(\omega)")
    ylabel("Im(\omega)")
    
    subplot(3,3,7); cla;
    h_spec_m1 = visualize_DMD_spectrum(out_METH1.omega, out_METH1.b, max(TimeDomain), "METH1"); hold on;
    h_spec_m2 = visualize_DMD_spectrum(out_METH2.omega, out_METH2.b, max(TimeDomain), "METH2"); hold off;
    [h_spec_m1.Color] = deal('r');
    [h_spec_m1.MarkerSize] = deal(8);
    
    [h_spec_m2.Color] = deal('b');
    [h_spec_m2.MarkerSize] = deal(8);
    legend('fontsize',6,'location','best');
    
    
    sgtitle({"Reduced-order models; TimeScale = " + TimeScale + " Duration = "+max(TimeDomain); ...
        IncludedComponents; ExcludedComponents;
        },'horizontalAlignment', 'left' );
    
    
    frame = getframe(h);
    writeVideo(v,frame);
    
end

close(v);
export_fig(filename, '-png',h);


% 初期化
clearvars;

% パラメーター
Lx = 4;% 計算領域
dLx = 0.01;% 要素サイズ
dt = 0.005;
nx = Lx / dLx + 1;% 要素数
c = 1;% 移流速度
t_max = 3;
iter_max = t_max / dt;
dt_pos = 0.02;% 画像出力間隔

% クーラン数
nu = c * dt / dLx;

% 座標の生成
x = 0 : dLx : Lx;

% 配列確保
q_exact = zeros(nx, 1);
q = zeros(nx, 1);
flux = zeros(nx - 1, 1);

% 初期条件の設定
for i =  1 : nx

    if abs(x(i) - 1) <= 10^-6

        q_exact(i) = 1;
        q(i) = 1;

    elseif x(i)  < 1

        q_exact(i) = 1;
        q(i) = 1;

    else

        q_exact(i) = 0;
        q(i) = 0;

    end

end

% 移流スキームの選択
method = '1st order UDS';% 1次精度風上差分法
% method = '2nd order CDS';% 2次精度中心差分法
% method = '2nd order UDS';% 2次精度風上差分法
% method = 'Lax-Friedrich';% Lax-Friedrich法
% method = 'Lax-Wendroff';% Lax-Wendroff法


% 時間発展
for iter = 1 : iter_max

    time = iter * dt;
    q_old = q;

    % 解析解
    for i = 1 : nx
        if abs(x(i) - (1 + c * time)) < 10^-6

            q_exact(i) = 1;

        elseif x(i) - (1 + c * time) < 0
            
            q_exact(i) = 1;
        
        else
        
            q_exact(i) = 0;
        
        end
    end

    % 数値解
    flux = CalcFlux(flux, q_old, method, nx, c, nu);% フラックスの計算

    for i = 2 : nx - 1

        q(i) = q_old(i) - (dt / dLx) * (flux(i) - flux(i - 1));

    end


    % コマンドウィンドウへの出力
    txt = ['iter = ', num2str(iter), ' / ', num2str(iter_max)];
    disp(txt);

    % リアルタイム可視化
    filename = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu),'.gif'];
    fignum = 1;
    plot(x, q_exact, x, q)
    title(['time = ', num2str(time, '%.3f')]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    axis equal; axis tight; axis on;
    fig=gcf;
    fig.Color='white';
    xlim([0 4]);
    ylim([-1 2]);
    xlabel('x')
    ylabel('q')

    legtxt = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu)];
    legend('exact', legtxt,'Location','southwest')
    legend('boxoff')

    frame = getframe(fignum);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if time == dt_pos
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.05, 'Loopcount', inf);
    elseif rem(time, dt_pos) == 0
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.05, 'WriteMode', 'append');
    end

end

function [flux] = CalcFlux(flux, q, method, nx, c, nu)

switch method

    %flux(i)はf_{i+1/2}に相当

    case '1st order UDS'

        for i = 1 : nx - 1

            flux(i) = c * q(i);% 1次精度風上差分

        end

    case '2nd order CDS'

        for i = 1 : nx - 1

            flux(i) = c * (q(i + 1) + q(i)) / 2;% 2次精度中心差分

        end

    case '2nd order UDS'

        for i = 2 : nx - 1

            flux(i) =  3 / 2 * c * q(i) - 1 / 2 * c * q(i - 1);% 2次精度風上差分

        end
        flux(1) = flux(2) ;

    case 'Lax-Friedrich'

        for i = 1 : nx - 1

            flux(i) = c / 2 * ((1 - 1 / nu) * q(i + 1) + (1 + 1 / nu) * q(i));% Lax-Friedrich法

        end

    case 'Lax-Wendroff'

        for i = 1 : nx - 1

            flux(i) = c / 2 * ((1 - nu) * q(i + 1) + (1 + nu) * q(i));% Lax-Wendroff法

        end

    otherwise

        error("選択された手法：" + method + "はサポートしていないです。" );
end

end



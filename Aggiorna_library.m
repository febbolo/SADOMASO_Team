%% PER AGGIORNARE LA LIBRERIA

% nuove versioni usando Source e non Reference, attento brodo
% basta cambiare old e new in base ai nomi che salvano

% se vuoi cambiare altro cambia anche il sistem da aprire

% Aggiorna la libreria per tutti i blocchi trovati

open_system('sim_lab7_disturbances.slx')

oldLib = 'LIBRARY';
newLib = 'LIBRARY_SAD';

blks = find_system(gcs, 'LookUnderMasks','all','FollowLinks','on');

for i = 1:numel(blks)
    try
        src = get_param(blks{i}, 'SourceBlock');

        % Deve iniziare per 'LIBRARY/'
        if startsWith(src, [oldLib '/'])
            % Rimpiazza SOLO il prefisso iniziale
            newSrc = [newLib src(length(oldLib)+1:end)];
            set_param(blks{i}, 'SourceBlock', newSrc);
        end

    catch
    end
end

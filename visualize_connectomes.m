%Este codigo permite visualizar todos los conectomas tanto del grupo control como del grupo experimental. 

control_folder = '/MATLAB Drive/Linux/Tractography/matriz_labels/control_connectomes';
filepatern = fullfile(control_folder, '*.csv'); 
files = dir(filepatern)

labels_table = readtable('Copy_of_label_fs_default.txt');
area_names = labels_table.Var2

%%%% Iniciar el ciclo for
for i = 1:length(files)
   current_filename = files(i).name;
   filepath = fullfile(control_folder, current_filename);

   connectome = importdata(filepath); 
   connectome_log = log1p(connectome);
   figure; 
   imagesc(connectome_log); 
   colorbar;
   colormap('parula'); 
   title(['Matriz de Conectividad: ' current_filename], 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
  
   xticks(1:length(area_names));
   xticklabels(area_names);
   xtickangle(90);
   yticks(1:length(area_names));
   yticklabels(area_names);
   set(gca, 'FontSize', 5);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experimental_folder = '/MATLAB Drive/Linux/Tractography/matriz_labels/experimental_connectomes';
filepaternschz = fullfile(schz_folder, '*.csv'); 
filesschz = dir(filepaternschz)

labels_table = readtable('/MATLAB Drive/Linux/Tractography/matriz_labels/Copy_of_label_fs_default.txt');
area_names = labels_table.Var2

%%%% Iniciar el ciclo for
for i = 1:length(filesschz)
    current_filename = filesschz(i).name;
    filepath = fullfile(schz_folder, current_filename);

    connectome = importdata(filepath); 
    connectome_log = log1p(connectome);
    figure; 
    imagesc(connectome_log); 
    colorbar;
    colormap('parula'); 
    %title([Matriz de Conectividad: ' current_filename], 'Interpreter', 'none', 'FontSize', 30);
    %title('\fontsize{12}Matriz de Conectividad' current_ ')

    xticks(1:length(area_names));
    xticklabels(area_names);
    xtickangle(90);
    yticks(1:length(area_names));
    yticklabels(area_names);
    set(gca, 'FontSize', 5);

    ti = title(['Matriz de Conectividad: ' current_filename], 'Interpreter', 'none', 'FontWeight', 'bold');
    ti.FontSize = 12;
end

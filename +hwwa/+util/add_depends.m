function add_depends()

conf = hwwa.config.load();

repos = conf.DEPENDS.repositories;
repo_p = conf.PATHS.repositories;
others = conf.DEPENDS.others;

for i = 1:numel(repos)
  addpath( genpath(fullfile(repo_p, repos{i})) );
end

for i = 1:numel(others)
  addpath( genpath(others{i}) );
end

end
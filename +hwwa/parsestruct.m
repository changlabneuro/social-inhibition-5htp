function params = parsestruct(params,args)

try
  params = shared_utils.general.parsestruct( params, args );
catch err
  throw( err );
end

end
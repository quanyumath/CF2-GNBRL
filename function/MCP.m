function y = MCP(x,alpha)

  if x <= alpha
      y = x-x.^2/2/alpha;
  else
      y = alpha/2;
  end
  
end

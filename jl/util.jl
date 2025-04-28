using GLMakie

function delete_all!(fig::Figure)
  while !isempty(fig.content)
    for child in fig.content delete!(fig.content[1]) end
    trim!(fig.layout)
  end
end


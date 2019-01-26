function print_block(fun::Function,out=stdout)
    io = IOBuffer()
    fun(io)
    output = String(take!(io))
    isempty(output) && return
    data = split(output, "\n")
    isempty(data[end]) && (data = data[1:end-1])
    if length(data) == 1
        println(out, "[ $(data[1])")
    elseif length(data) > 1
        for (p,dl) in zip(vcat("╭", repeat(["│"], length(data)-2), "╰"),data)
            println(out, "$(p) $(dl)")
        end
    end
    println(out)
end

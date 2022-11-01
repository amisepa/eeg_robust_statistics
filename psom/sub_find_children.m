% GRAPH_DEPS(J,K) == 1 if and only if JOB K depends on JOB J. GRAPH_DEPS =
% 0 otherwise. This (ugly but reasonably fast) recursive code will work
% only if the directed graph defined by GRAPH_DEPS is acyclic.
% MASK_CHILD(NUM_J) == 1 if the job NUM_J is a children of one of the job
% in the boolean mask MASK and the job is in MASK_TODO.
% This last restriction is used to speed up computation.

function mask_child = sub_find_children(mask,graph_deps)

if max(double(mask))>0
    mask_child = max(graph_deps(mask,:),[],1)>0;    
    mask_child_strict = mask_child & ~mask;
else
    mask_child = false(size(mask));
end
if any(mask_child)
    mask_child = mask_child | sub_find_children(mask_child_strict,graph_deps);
end

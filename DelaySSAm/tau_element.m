function result = tau_element(element)
    % TAU_ELEMENT Checks if the input is a function handle.
    % If it is a function handle, it evaluates the function.
    % Otherwise, it returns the input element as is.
    %
    % Syntax:
    % result = tau_element(element)
    %
    % Inputs:
    % element - A function handle or any other type of input.
    %
    % Outputs:
    % result - The result of the function evaluation or the input element itself.

    if isa(element, 'function_handle')
        result = element();
    else
        result = element;
    end
end
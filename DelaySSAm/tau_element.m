function result = tau_element(element)

    if isa(element, 'function_handle')
        result = element();
    else
        result = element;
    end
end
"""
# Since in aiida-fleur we are parsing and workchains are reacting on several errors 
# reported by fleur which we recognize by the 'error message', this interface has to be tested and
# should stay the same if possible (otherwise this might lead to back compability issues 
# to previous fleur versions).
# Or at least if the messages are refactored, i.e one of these tests here fails
# aiida-fleur has to be maintained. In this case please inform aiida-fleur devs
"""
import pytest


def test_judft_messages_aiida_interface(collect_all_judft_messages_f):
    """
    While the same judft messages can be thrown in different places in Fleur, 
    we test if the 'message identifier' which is used inside aiida-fleur 
    if expected uniquely identifies with a judft case (message).
    We test only messages critical for aiida-fleur here, but maybe one wants to extent this.

    One could also regression test the whole judft_message interface if wanted.
    """
    
    # TODO import this from somewhere. (masci-tools?)
    # tuple, (messagestring, expected_unique)
    aiida_fleur_parsed_messages = [
        #(r'cgroup out-of-memory handler', True), # This is not in fleur
        #('Out Of Memory', True),  # This is not in fleur
        ('Allocation of array for communication failed', True),
        ('Atom spills out into vacuum during relaxation', True),
        ('Error checking M.T. radii', True),
        ('Overlapping MT-spheres during relaxation: ', True),
        ('FLEUR is not linked against libxc', True),
        ('No solver linked for Hubbard 1', True)]

    judft_error_message_list = collect_all_judft_messages_f()
    print(judft_error_message_list)
    for message, unique in aiida_fleur_parsed_messages:
        print(message)
        count = -1
        dev_msg = (f'The judft message part: "{message}" has changed,'
                   ' if this interface change is necessary, please inform an aiida-fleur dev.')
        for judft_msg in judft_error_message_list:
            if message in str(judft_msg):
                count = count +1
        if count >0 and unique:
            dev_msg = dev_msg + " The message is not uniquely identifiable anymore :'(. This maybe wanted but please, check."
            assert False, dev_msg
        if count == -1:
            dev_msg = dev_msg + " The message identifier does not exists anymore :'(."
            assert False, dev_msg

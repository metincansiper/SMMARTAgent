from tumor_board_agent.tumor_board_agent import TumorBoardAgent

def join_list(l):
    return ','.join(l)

tba = TumorBoardAgent()
while True:
    print( 'Create a tumor board by:' )
    print( '1) Patient id' )
    print( '2) List of variant genes' )
    print( '3) Exit' )
    param_type = input(' >> ')

    prompt_param_str = None

    if param_type == '1':
        prompt_param_str = '\nPlease enter the patient id : '
    elif param_type == '2':
        prompt_param_str = '\nPlease enter the list of variant genes seperated by space : '
    elif param_type == '3':
        break
    else:
        # raise Exception('Invalid input {}'.format(param_type))
        print('\nInvalid input {}'.format(param_type))
        continue

    param = input(prompt_param_str)

    if param_type == '2':
        param = param.split()

    results = tba.create_tumor_board_report(param)
    print('\nI found {} genes in the tumor board report : '.format(len(results)))
    print('\n')
    print(join_list(results))
    # print(results)

    while True:
        print('\n')
        print('1) List the top genes')
        print('2) Ask why a gene is interesting')
        print('3) Ask for evidence for interaction')
        print('4) Go back to main menu')

        choice = input(' >> ')

        if choice == '1':
            k = input('\nPlease enter number of top genes to list : ')
            k = int(k)
            top_genes = tba.get_top_k_res(k)
            print(join_list(top_genes))
        elif choice == '2':
            gene = input('\nPlease enter the name of gene of interest : ')
            why_important = tba.why_important(gene, None)
            print('\n{} is interesting because it interacts with mutated genes including {}'.format(gene, join_list(why_important)))
        elif choice == '3':
            gene1 = input('\nPlease enter a gene of the interaction : ')
            gene2 = input('\nPlease enter the other gene of the interaction : ')

            evidences = tba.get_evidences_for(gene1, gene2)
            SHORT_LIST_SIZE = 5
            short_list = evidences[:SHORT_LIST_SIZE]

            print('\n{} interacts with {} because Pathway Commons has it.'.format(gene1, gene2))
            print('\nI found following evidences for this interaction: {}'.format(join_list(short_list)))

            if len(evidences) > SHORT_LIST_SIZE:
                print('\nI found more evidences. Would you like to view the whole list?')
                print('1) Yes')
                print('2) No')
                sub_choice = input(' >> ')

                if sub_choice == '1':
                    print(join_list(evidences))

        elif choice == '4':
            break
        else:
            # raise Exception('Invalid input {}'.format(choice))
            print('\nInvalid input {}'.format(choice))
            continue

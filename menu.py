from tumor_board_agent.tumor_board_agent import TumorBoardAgent

tba = TumorBoardAgent()
while True:
    print( 'Create tumor board by:' )
    print( '1) Patient id' )
    print( '2) List of variant genes' )
    print( '3) Exit' )
    param_type = input(' >> ')

    prompt_param_str = None

    if param_type == '1':
        prompt_param_str = 'Please enter patient id : '
    elif param_type == '2':
        prompt_param_str = 'Please enter list of variant genes seperated by space : '
    elif param_type == '3':
        break
    else:
        # raise Exception('Invalid input {}'.format(param_type))
        print('Invalid input {}'.format(param_type))
        continue

    param = input(prompt_param_str)

    if param_type == '2':
        param = param.split()

    results = tba.create_tumor_board_report(param)
    print(results)

    while True:
        print('1) List top k genes')
        print('2) Ask why a gene is important')
        print('3) Ask for evidence for interaction')
        print('4) Go back to main menu')

        choice = input(' >> ')

        if choice == '1':
            k = input('Please enter k : ')
            k = int(k)
            top_genes = tba.get_top_k_res(k)
            print(top_genes)
        elif choice == '2':
            gene = input('Please enter the name of gene of interest : ')
            why_important = tba.why_important(gene)
            print('{} is important because it interacts with mutated genes including {}'.format(gene, ','.join(why_important)))
        elif choice == '3':
            gene1 = input('Please enter one gene of the interaction : ')
            gene2 = input('Please enter the other gene of the interaction : ')

            evidences = tba.get_evidences_for(gene1, gene2)
            SHORT_LIST_SIZE = 5
            short_list = evidences[:SHORT_LIST_SIZE]

            print('{} interacts with {} because Pathway Commons has it.'.format(gene1, gene2))
            print('We have the fallowing evidences for this interaction: {}'.format(','.join(short_list)))

            if len(evidences) > SHORT_LIST_SIZE:
                print('We have found more evidences. Would you like to view the whole list?')
                print('1) Yes')
                print('2) No')
                sub_choice = input(' >> ')

                if sub_choice == '1':
                    print(','.join(evidences))

        elif choice == '4':
            break
        else:
            # raise Exception('Invalid input {}'.format(choice))
            print('Invalid input {}'.format(choice))
            continue

from tumor_board_agent.tumor_board_agent import TumorBoardAgent

sbgn = TumorBoardAgent.get_pathsbetween_genes(['TP53', 'MDM2'])
print(sbgn)

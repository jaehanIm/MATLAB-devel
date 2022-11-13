function setImuIdx(jobIdx)
    global missStartTime missEndTime start finish imuTimeS
    gdStart = missStartTime(jobIdx);
    gdFinish = missEndTime(jobIdx);
    start = min(find(imuTimeS>gdStart & imuTimeS<gdFinish));
    finish = max(find(imuTimeS>gdStart & imuTimeS<gdFinish));
end
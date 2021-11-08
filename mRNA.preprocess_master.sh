#!/bin/bash
# @author Slim Fourati (slim.fourati@emory.edu)
# @version 0.8

# read input arguments
email="slim.fourati@emory.edu"
genome=GRCh38
pairEnd=false
isoform=false
acceptedGenome=("GRCh38" "Mmul_10" "Mnem_1")

while getopts :d:e:g:pih option
do
    case "${option}" in
	h) echo "Command: bash mRNA.preprocess_master.sh -d {fastq/directoryfastq} ..."
	    echo "argument: d=[d]irectory with raw data (required)"
	    echo "          g=reference [g]enome"
	    echo "          p=[p]aired-end sequencing"
	    echo "          i=[i]soform transcript/exon counts"
	    echo "          h=print [h]elp"
	    exit 1;;
	d) dirFastq=$OPTARG;;
	e) email=$OPTARG;;
	g) genome=$OPTARG
	    if [[ ! "${acceptedGenome[@]}" =~ "$genome" ]]
	    then
		echo "Invalid -g argument: choose between ${acceptedGenome[@]}"
		exit 1
	    fi;;
	p) pairEnd=true;;
	i) isoform=true;;
	\?) echo "Invalid option: -$OPTARG"
	    exit 1;;
	:)
	    echo "Option -$OPTARG requires an argument."
	    exit 1;;
    esac
done

# test that directory is provided
if [ -z ${dirFastq+x} ]
then
    echo "error...option -d required."
    exit 1
fi

# test that directory contains seq files
suffix="fq.gz"
nfiles=$(find $dirFastq -name "*_1.$suffix" | wc -l)
if [ $nfiles -lt 1 ]
then
    echo "error...empty input directory"
    exit 1
fi

# initialize directory
dirData=$(echo $dirFastq | sed -r "s|efs/||g")

# lauch genome indexing
sed -ri "s|\"jobName\": \"rnaseq-job-TIMESTAMP\",|\"jobName\": \"rnaseq-job-$(date +%Y%m%d%M%S)\",|g" \
    genomeGenerate.json
sed -ri "s|\"-d\",.+$|\"-d\",\"${dirData}\",|g" \
    genomeGenerate.json
sed -ri "s|\"-g\",.+$|\"-g\",\"${genome}\"|g" \
    genomeGenerate.json

# copy script on mount
cp genomeGenerate.sh /mnt/efs/

cmd="aws batch submit-job"
cmd="$cmd --cli-input-json file://genomeGenerate.json"
cmd="$cmd --profile 'tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole'"

# echo $cmd
jobid=$(eval $cmd | \
	    grep jobId | \
	    sed -r 's|.+jobId\": \"(.+)\"$|\1|g')

# modify rule.json
ruleName="rule-$(date +%Y%m%d%M%S)"
sed -ri "s|\"Name\":.+$|\"Name\": \"${ruleName}\",|g" \
    rule.json
sed -ri "s|\"jobId\\\\\":[^,]+|\"jobId\\\\\":[\\\\\"${jobid}\\\\\"]|g" \
    rule.json
sed -ri "s|\"Description\":.+$|\"Description\": \"${jobid}\",|g" \
    rule.json

# register rule
aws events put-rule \
    --cli-input-json file://rule.json \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# register topic
cmd="aws sns create-topic"
cmd="$cmd --name 'job-$(date +%Y%m%d%M%S)'"
cmd="$cmd --profile 'tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole'"
topicArn=$(eval $cmd | \
            grep TopicArn | \
            sed -r 's|.+TopicArn\": \"(.+)\"$|\1|g')

# add iam role to allow notification
attributeValue="{\"Version\":\"2012-10-17\",\"Id\":\"__default_policy_ID\",\"Statement\":[{\"Sid\":\"__default_statement_ID\",\"Effect\":\"Allow\",\"Principal\":{\"AWS\":\"*\"},\"Action\":[\"SNS:GetTopicAttributes\",\"SNS:SetTopicAttributes\",\"SNS:AddPermission\",\"SNS:RemovePermission\",\"SNS:DeleteTopic\",\"SNS:Subscribe\",\"SNS:ListSubscriptionsByTopic\",\"SNS:Publish\",\"SNS:Receive\"],\"Resource\":\"${topicArn}\",\"Condition\":{\"StringEquals\":{\"AWS:SourceOwner\":\"943708588556\"}}},{\"Sid\":\"AWSEvents_${ruleName}_1\",\"Effect\":\"Allow\",\"Principal\":{\"Service\":\"events.amazonaws.com\"},\"Action\":\"sns:Publish\",\"Resource\":\"${topicArn}\"}]}"
aws sns set-topic-attributes \
    --topic-arn "$topicArn" \
    --attribute-name "Policy" \
    --attribute-value "$attributeValue" \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# add target to rule
aws events put-targets \
    --rule "$ruleName" \
    --targets "Id"=1,"Arn"="$topicArn" \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# create subscription
aws sns subscribe \
    --topic-arn "$topicArn" \
    --protocol "email" \
    --notification-endpoint "$email" \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# modify preprocessing json
sed -ri "s|\"jobName\": \"rnaseq-job-TIMESTAMP\",|\"jobName\": \"rnaseq-job-$(date +%Y%m%d%M%S)\",|g" \
    mRNA.preprocess_seq.json
sed -ri "s|\"jobId\":.+$|\"jobId\": \"${jobid}\"|g" \
    mRNA.preprocess_seq.json
sed -ri "s|\"size\":.+$|\"size\": ${nfiles}|g" \
    mRNA.preprocess_seq.json
sed -ri "s|\"-d\",.+$|\"-d\",\"${dirData}\",|g" \
    mRNA.preprocess_seq.json
sed -ri "s|\"-g\",.+$|\"-g\",\"${genome}\",|g" \
    mRNA.preprocess_seq.json
if ! $pairEnd
then
    sed -rzi "s|,\n[^\n]+\"-p\"||g" mRNA.preprocess_seq.json
fi
if ! $isoform
then
    sed -rzi "s|,\n[^\n]+\"-i\"||g" mRNA.preprocess_seq.json
fi

# lauch preprocessing script
cmd="aws batch submit-job"
cmd="$cmd --cli-input-json file://mRNA.preprocess_seq.json"
cmd="$cmd --profile 'tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole'"

# copy script on mount
cp mRNA.preprocess_seq.sh /mnt/efs/

# echo $cmd
jobid=$(eval $cmd | \
            grep jobId | \
            sed -r 's|.+jobId\": \"(.+)\"$|\1|g')

# modify rule.json
ruleName="rule-$(date +%Y%m%d%M%S)"
sed -ri "s|\"Name\":.+$|\"Name\": \"${ruleName}\",|g" \
    rule.json
sed -ri "s|\"jobId\\\\\":[^,]+|\"jobId\\\\\":[\\\\\"${jobid}\\\\\"]|g" \
    rule.json
sed -ri "s|\"Description\":.+$|\"Description\": \"${jobid}\",|g" \
    rule.json

# register rule
aws events put-rule \
    --cli-input-json file://rule.json \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# add iam role to allow notification
attributeValue=$(echo $attributeValue |
		     sed -r "s|]}$||g")
attributeValue="${attributeValue},{\"Sid\":\"AWSEvents_${ruleName}_1\",\"Effect\":\"Allow\",\"Principal\":{\"Service\":\"events.amazonaws.com\"},\"Action\":\"sns:Publish\",\"Resource\":\"${topicArn}\"}]}"
aws sns set-topic-attributes \
    --topic-arn "$topicArn" \
    --attribute-name "Policy" \
    --attribute-value "$attributeValue" \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"

# add target to rule
aws events put-targets \
    --rule "$ruleName" \
    --targets "Id"=1,"Arn"="$topicArn" \
    --profile "tki-aws-account-310-rhedcloud/RHEDcloudAdministratorRole"
